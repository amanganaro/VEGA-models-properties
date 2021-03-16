package insilico.carcinogenicity_sfiregression.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.HashMap;
import java.util.Map;

/**
 * Weighting scheme for EState. Calculation taken from code by EPA (T.E.S.T)
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class EState {


    private Map<String, Integer> periods;

    private final double[] D;
    private final double[] DV;
    private final double[] DV2;
    private final double[] IS;
    private final double[] EState;
    private final double[] HEState;
    private final double[] KHE;

    public EState(IAtomContainer mol) throws GenericFailureException {

        SetPeriods();

        int nSk = mol.getAtomCount();

        // Get matrices
        int[][] TopDistMat;
        try {
            TopDistMat = TopoDistanceMatrix.getMatrix(mol);
        } catch (Exception e) {
            log.warn(e.getMessage());
            throw new GenericFailureException("Unable to calculate matrices");
        }

        // Initialize variables
        D = new double[nSk];
        DV = new double[nSk];
        DV2 = new double[nSk];
        IS = new double[nSk];
        EState = new double[nSk];
        HEState = new double[nSk];
        KHE = new double[nSk];

        for (int i = 0; i<nSk; i++) {
            D[i] = 0.0;
            DV[i] = 0;
            DV2[i] = 0.0;
            IS[i] = 0.0;
            EState[i] = 0;
            HEState[i] = 0;
            KHE[i] = 0;
        }

        // Calculations
        CalculateDeltaValues(mol);
        CalculateDV(mol);
        CalculateIntrinsicStateAndKHE(mol);
        CalculateEState(mol, TopDistMat);
    }

    private void CalculateDeltaValues(IAtomContainer mol) throws GenericFailureException {

        // determines D, DV2 for each specific atom
        for (int I = 0; I <= mol.getAtomCount() - 1; I++) {

            //DV2[I]=m.getBondOrderSum(m.getAtom(I));
            //modified to take into account aromaticity
            //to avoid problems with smiles that give
            //different values (es c1ccccc1 - C1=CC=CC1)
            double count = 0;
            for (IBond bond : mol.getConnectedBondsList(mol.getAtom(I))) {
                if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                    count += 1.5;
                } else {
                    if (bond.getOrder() == IBond.Order.SINGLE) {
                        count += 1.0;
                    } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                        count += 2.0;
                    } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                        count += 3.0;
                    } else if (bond.getOrder() == IBond.Order.QUADRUPLE) {
                        count += 4.0;
                    }
                }
            }

            DV2[I]=count;
            D[I]=mol.getConnectedBondsCount(mol.getAtom(I));

            String symbol = mol.getAtom(I).getSymbol();
            int charge = mol.getAtom(I).getFormalCharge();

            if (symbol.equals("C") || symbol.equals("Si")
                    || symbol.equals("Pb") || symbol.equals("Sn"))
            {

                // TODO: I am not sure if following block is needed:
                if ((int) (getD()[I]) != (int) (0.5 + getDV2()[I])) {

                    //used for carbons with fractional aromatic bond order
                    int ID = (int) getD()[I];

                    switch (ID) {
                        case 1: { // 41
                            if ((int) (getDV2()[I] + 0.5) != 2)
                                DV2[I] = 3.0;
                            break;
                        }
                        case 2: { // 42
                            if ((int) (getDV2()[I] + 0.5) != 3)
                                DV2[I] = 4.0;
                            break;
                        }
                        case 3: { // 43
                            if (getDV2()[I] >= 4)
                                DV2[I] = 4.0;
                            break;
                        }
                    }
                }

            } else if (symbol.equals("O")) { // 50 (oxygen)

                if (getDV2()[I] > 1)
                    DV2[I] = 6.0;
                else
                    DV2[I] = 5.0;

            } else if (symbol.equals("N")) { // 60 (nitrogen)

                int ID = (int) getD()[I];
                switch (ID) {
                    case 1: // 61-sp3
                    {
                        DV2[I] = getDV2()[I] + 2.;
                        break;
                    }
                    case 2: // 62-sp2
                    {
                        if ((int) (getDV2()[I] + 0.5) == 2)
                            DV2[I] = 4.;
                        else
                            DV2[I] = 5.;
                        break;
                    }
                    case 3: // 63-sp
                    {
                        if ((int) (getDV2()[I] + 0.5) == 3 || (int) (getDV2()[I] + 0.5) == 4)
                            DV2[I] = getDV2()[I] + 2.;
                        else
                            DV2[I] = 5.;
                        break;
                    }
                }

            } else if (symbol.equals("S")) { // 70 (sulfur)
                // ......SULFUR

                double h,Zv;
                Zv=6;
                if (getDV2()[I]==1)
                    h=1;
                else
                    h=0;
                DV2[I]=Zv-h;

            } else if (symbol.equals("P") || symbol.equals("As")) {

                double h,Zv;
                Zv=5;

                if (getDV2()[I]==1)
                    h=2;
                else if (getDV2()[I]==2)
                    h=1;
                else if (getDV2()[I]==4)
                    h=1;
                else
                    h=0;
                DV2[I]=Zv-h;

            } else { // F,Cl,Br,I, Hg

                double h = 0; // assume number of hydrogens = 0
                double Zv = ValenceVertexDegree.GetValenceElectronsNumber(symbol);
                if (Zv == Descriptor.MISSING_VALUE)
                    throw new GenericFailureException("unable to set Zv for atom type: " + symbol);

                DV2[I] = (Zv - h);

            }

            // modify based on charge
            DV2[I] += -charge; // alternatively could use the charge to modify the number of hydrogens
        }
    }

    private void CalculateDV(IAtomContainer mol) throws GenericFailureException {

        for (int i = 0; i <= mol.getAtomCount() - 1; i++) {
            String symbol=mol.getAtom(i).getSymbol();

            double Zv = ValenceVertexDegree.GetValenceElectronsNumber(symbol);
            if (Zv == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("unable to set Zv for atom type: " + symbol);
            double Z = (double)mol.getAtom(i).getAtomicNumber();

            DV[i]=getDV2()[i]/(Z - Zv - 1);
        }
    }


    private void CalculateIntrinsicStateAndKHE(IAtomContainer mol) {
        double Z=0, Zv=0;
        for (int i = 0; i <= mol.getAtomCount() - 1; i++) {
            IAtom a = mol.getAtom(i);

            int period = periods.get(a.getSymbol());
            double N = (double)period;


            IS[i] = (Math.pow(2.0 / N, 2.0) * getDV2()[i] + 1) / getD()[i];
            KHE[i] = (getDV2()[i] - getD()[i]) / (N * N);
        }
    }


    private void CalculateEState(IAtomContainer mol, int[][] DistanceMatrix) {

        for (int i = 0; i <= mol.getAtomCount() - 1; i++) {
            double sumDeltaIijForE = 0;
            double sumDeltaIijForHE = 0;

            for (int j = 0; j <= mol.getAtomCount() - 1; j++)
                if (i != j) {
                    sumDeltaIijForE += (getIS()[i] - getIS()[j])/ Math.pow((double) DistanceMatrix[i][j] + 1.0, 2.0);
                    sumDeltaIijForHE += (getKHE()[j] + 0.2)/ Math.pow((double) DistanceMatrix[i][j] + 1.0, 2.0);
                }

            EState[i] = getIS()[i] + sumDeltaIijForE;
            HEState[i] = getKHE()[i] + (getKHE()[i] + 0.2) + sumDeltaIijForHE;
        }
    }


    private void SetPeriods() {
        periods = new HashMap<>();
        periods.put("H", 1);
        periods.put("He", 1);
        periods.put("Li", 2);
        periods.put("Be", 2);
        periods.put("B", 2);
        periods.put("C", 2);
        periods.put("N", 2);
        periods.put("O", 2);
        periods.put("F", 2);
        periods.put("Ne", 2);
        periods.put("Na", 3);
        periods.put("Mg", 3);
        periods.put("Al", 3);
        periods.put("Si", 3);
        periods.put("P", 3);
        periods.put("S", 3);
        periods.put("Cl", 3);
        periods.put("Ar", 3);
        periods.put("K", 4);
        periods.put("Ca", 4);
        periods.put("Sc", 4);
        periods.put("Ti", 4);
        periods.put("V", 4);
        periods.put("Cr", 4);
        periods.put("Mn", 4);
        periods.put("Fe", 4);
        periods.put("Co", 4);
        periods.put("Ni", 4);
        periods.put("Cu", 4);
        periods.put("Zn", 4);
        periods.put("Ga", 4);
        periods.put("Ge", 4);
        periods.put("As", 4);
        periods.put("Se", 4);
        periods.put("Br", 4);
        periods.put("Kr", 4);
        periods.put("Rb", 5);
        periods.put("Sr", 5);
        periods.put("Y", 5);
        periods.put("Zr", 5);
        periods.put("Nb", 5);
        periods.put("Mo", 5);
        periods.put("Tc", 5);
        periods.put("Ru", 5);
        periods.put("Rh", 5);
        periods.put("Pd", 5);
        periods.put("Ag", 5);
        periods.put("Cd", 5);
        periods.put("In", 5);
        periods.put("Sn", 5);
        periods.put("Sb", 5);
        periods.put("Te", 5);
        periods.put("I", 5);
        periods.put("Xe", 5);
        periods.put("Cs", 6);
        periods.put("Ba", 6);
        periods.put("La", 6);
        periods.put("Ce", 6);
        periods.put("Pr", 6);
        periods.put("Nd", 6);
        periods.put("Pm", 6);
        periods.put("Sm", 6);
        periods.put("Eu", 6);
        periods.put("Gd", 6);
        periods.put("Tb", 6);
        periods.put("Dy", 6);
        periods.put("Ho", 6);
        periods.put("Er", 6);
        periods.put("Tm", 6);
        periods.put("Yb", 6);
        periods.put("Lu", 6);
        periods.put("Hf", 6);
        periods.put("Ta", 6);
        periods.put("W", 6);
        periods.put("Re", 6);
        periods.put("Os", 6);
        periods.put("Ir", 6);
        periods.put("Pt", 6);
        periods.put("Au", 6);
        periods.put("Hg", 6);
        periods.put("Tl", 6);
        periods.put("Pb", 6);
        periods.put("Bi", 6);
        periods.put("Po", 6);
        periods.put("At", 6);
        periods.put("Rn", 6);
        periods.put("Fr", 7);
        periods.put("Ra", 7);
        periods.put("Ac", 7);
        periods.put("Th", 7);
        periods.put("Pa", 7);
        periods.put("U", 7);
        periods.put("Np", 7);
        periods.put("Pu", 7);
        periods.put("Am", 7);
        periods.put("Cm", 7);
        periods.put("Bk", 7);
        periods.put("Cf", 7);
        periods.put("Es", 7);
        periods.put("Fm", 7);
        periods.put("Md", 7);
        periods.put("No", 7);
        periods.put("Lr", 7);
        periods.put("Rf", 7);
        periods.put("Db", 7);
        periods.put("Sg", 7);
        periods.put("Bh", 7);
        periods.put("Hs", 7);
        periods.put("Mt", 7);
        periods.put("Ds", 7);
        periods.put("Rg", 7);
        periods.put("Uub", 7);
        periods.put("Uut", 7);
        periods.put("Uuq", 7);
        periods.put("Uup", 7);
        periods.put("Uuh", 7);
        periods.put("Uus", 7);
        periods.put("Uuo", 7);
    }

    /**
     * @return the D
     */
    public double[] getD() {
        return D;
    }

    /**
     * @return the DV
     */
    public double[] getDV() {
        return DV;
    }

    /**
     * @return the DV2
     */
    public double[] getDV2() {
        return DV2;
    }

    /**
     * @return the IS
     */
    public double[] getIS() {
        return IS;
    }

    /**
     * @return the EState
     */
    public double[] getEState() {
        return EState;
    }

    /**
     * @return the HEState
     */
    public double[] getHEState() {
        return HEState;
    }

    /**
     * @return the KHE
     */
    public double[] getKHE() {
        return KHE;
    }


}
