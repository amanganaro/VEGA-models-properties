//
// correzione sul calcolo di DV2 e D per la presenza di H
// quando viene chiamato da desc H-filled come autocorrelation
//

package insilico.aromatase_irfmn.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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

public class EStateCorrectForH {

    private static final Logger log = LogManager.getLogger(EStateCorrectForH.class);

    private Map<String, Integer> periods;

    private final double [] D;
    private final double [] DV;
    private final double [] DV2;
    private final double []IS;
    private final double []EState;
    private final double []HEState;
    private final double[] KHE;


    public EStateCorrectForH(IAtomContainer m) throws GenericFailureException {

        SetPeriods();

        // Get structure
//        Molecule m;
//        try {
//            m = Mol.GetStructure();
//        } catch (InvalidMoleculeException e) {
//            throw new GenericFailureException("Invalid molecule");
//        }
        int nSk = m.getAtomCount();

        // Get matrices
        int[][] TopoDistMat;
        try {
//            TopoDistMat = Mol.GetMatrixTopologicalDistance();
            TopoDistMat = TopoDistanceMatrix.getMatrix(m);
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
        CalculateDeltaValues(m);
        CalculateDV(m);
        CalculateIntrinsicStateAndKHE(m);
        CalculateEState(m, TopoDistMat);

    }


    private void CalculateDeltaValues(IAtomContainer m) throws GenericFailureException {
        // determines D, DV2 for each specific atom

        for (int I = 0; I <= m.getAtomCount() - 1; I++) {

            int nH = 0;
            for (IAtom a : m.getConnectedAtomsList(m.getAtom(I))) {
                if (a.getAtomicNumber() == 1)
                    nH++;
            }

            //DV2[I]=m.getBondOrderSum(m.getAtom(I));
            //modified to take into account aromaticity
            //to avoid problems with smiles that give
            //different values (es c1ccccc1 - C1=CC=CC1)
            double count = 0;
            for (IBond bond : m.getConnectedBondsList(m.getAtom(I))) {
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

            DV2[I]=count - nH;
            D[I]=m.getConnectedBondsCount(m.getAtom(I)) - nH;

            String symbol = m.getAtom(I).getSymbol();
            int charge = m.getAtom(I).getFormalCharge();

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


    private void CalculateDV(IAtomContainer m) throws GenericFailureException {

        for (int i = 0; i <= m.getAtomCount() - 1; i++) {
            String symbol=m.getAtom(i).getSymbol();

            double Zv = ValenceVertexDegree.GetValenceElectronsNumber(symbol);
            if (Zv == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("unable to set Zv for atom type: " + symbol);
            double Z = (double)m.getAtom(i).getAtomicNumber();

            DV[i]=getDV2()[i]/(Z - Zv - 1);
        }
    }


    private void CalculateIntrinsicStateAndKHE(IAtomContainer m) {
        double Z=0, Zv=0;
        for (int i = 0; i <= m.getAtomCount() - 1; i++) {
            IAtom a = (IAtom)m.getAtom(i);

            int period=(Integer)(periods.get(a.getSymbol()));
            double N = (double)period;

//            double VertexDegree = 0;
//            int H = Manipulator.CountImplicitHydrogens(a);
//            for (IAtom connAt : m.getConnectedAtomsList(a))
//                if (!connAt.getSymbol().equalsIgnoreCase("H"))
//                    VertexDegree++;
//            double ValenceVertexDegree = VertexDegree - H - a.getFormalCharge();
//            IS[i] = (Math.pow(2.0 / N, 2.0) * ValenceVertexDegree + 1) / VertexDegree;

            IS[i] = (Math.pow(2.0 / N, 2.0) * getDV2()[i] + 1) / getD()[i];
            KHE[i] = (getDV2()[i] - getD()[i]) / (N * N);
        }
    }


    private void CalculateEState(IAtomContainer m, int[][] DistanceMatrix) {

        for (int i = 0; i <= m.getAtomCount() - 1; i++) {
            double sumDeltaIijForE = 0;
            double sumDeltaIijForHE = 0;

            for (int j = 0; j <= m.getAtomCount() - 1; j++)
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
        periods.put("H", new Integer(1));
        periods.put("He", new Integer(1));
        periods.put("Li", new Integer(2));
        periods.put("Be", new Integer(2));
        periods.put("B", new Integer(2));
        periods.put("C", new Integer(2));
        periods.put("N", new Integer(2));
        periods.put("O", new Integer(2));
        periods.put("F", new Integer(2));
        periods.put("Ne", new Integer(2));
        periods.put("Na", new Integer(3));
        periods.put("Mg", new Integer(3));
        periods.put("Al", new Integer(3));
        periods.put("Si", new Integer(3));
        periods.put("P", new Integer(3));
        periods.put("S", new Integer(3));
        periods.put("Cl", new Integer(3));
        periods.put("Ar", new Integer(3));
        periods.put("K", new Integer(4));
        periods.put("Ca", new Integer(4));
        periods.put("Sc", new Integer(4));
        periods.put("Ti", new Integer(4));
        periods.put("V", new Integer(4));
        periods.put("Cr", new Integer(4));
        periods.put("Mn", new Integer(4));
        periods.put("Fe", new Integer(4));
        periods.put("Co", new Integer(4));
        periods.put("Ni", new Integer(4));
        periods.put("Cu", new Integer(4));
        periods.put("Zn", new Integer(4));
        periods.put("Ga", new Integer(4));
        periods.put("Ge", new Integer(4));
        periods.put("As", new Integer(4));
        periods.put("Se", new Integer(4));
        periods.put("Br", new Integer(4));
        periods.put("Kr", new Integer(4));
        periods.put("Rb", new Integer(5));
        periods.put("Sr", new Integer(5));
        periods.put("Y", new Integer(5));
        periods.put("Zr", new Integer(5));
        periods.put("Nb", new Integer(5));
        periods.put("Mo", new Integer(5));
        periods.put("Tc", new Integer(5));
        periods.put("Ru", new Integer(5));
        periods.put("Rh", new Integer(5));
        periods.put("Pd", new Integer(5));
        periods.put("Ag", new Integer(5));
        periods.put("Cd", new Integer(5));
        periods.put("In", new Integer(5));
        periods.put("Sn", new Integer(5));
        periods.put("Sb", new Integer(5));
        periods.put("Te", new Integer(5));
        periods.put("I", new Integer(5));
        periods.put("Xe", new Integer(5));
        periods.put("Cs", new Integer(6));
        periods.put("Ba", new Integer(6));
        periods.put("La", new Integer(6));
        periods.put("Ce", new Integer(6));
        periods.put("Pr", new Integer(6));
        periods.put("Nd", new Integer(6));
        periods.put("Pm", new Integer(6));
        periods.put("Sm", new Integer(6));
        periods.put("Eu", new Integer(6));
        periods.put("Gd", new Integer(6));
        periods.put("Tb", new Integer(6));
        periods.put("Dy", new Integer(6));
        periods.put("Ho", new Integer(6));
        periods.put("Er", new Integer(6));
        periods.put("Tm", new Integer(6));
        periods.put("Yb", new Integer(6));
        periods.put("Lu", new Integer(6));
        periods.put("Hf", new Integer(6));
        periods.put("Ta", new Integer(6));
        periods.put("W", new Integer(6));
        periods.put("Re", new Integer(6));
        periods.put("Os", new Integer(6));
        periods.put("Ir", new Integer(6));
        periods.put("Pt", new Integer(6));
        periods.put("Au", new Integer(6));
        periods.put("Hg", new Integer(6));
        periods.put("Tl", new Integer(6));
        periods.put("Pb", new Integer(6));
        periods.put("Bi", new Integer(6));
        periods.put("Po", new Integer(6));
        periods.put("At", new Integer(6));
        periods.put("Rn", new Integer(6));
        periods.put("Fr", new Integer(7));
        periods.put("Ra", new Integer(7));
        periods.put("Ac", new Integer(7));
        periods.put("Th", new Integer(7));
        periods.put("Pa", new Integer(7));
        periods.put("U", new Integer(7));
        periods.put("Np", new Integer(7));
        periods.put("Pu", new Integer(7));
        periods.put("Am", new Integer(7));
        periods.put("Cm", new Integer(7));
        periods.put("Bk", new Integer(7));
        periods.put("Cf", new Integer(7));
        periods.put("Es", new Integer(7));
        periods.put("Fm", new Integer(7));
        periods.put("Md", new Integer(7));
        periods.put("No", new Integer(7));
        periods.put("Lr", new Integer(7));
        periods.put("Rf", new Integer(7));
        periods.put("Db", new Integer(7));
        periods.put("Sg", new Integer(7));
        periods.put("Bh", new Integer(7));
        periods.put("Hs", new Integer(7));
        periods.put("Mt", new Integer(7));
        periods.put("Ds", new Integer(7));
        periods.put("Rg", new Integer(7));
        periods.put("Uub", new Integer(7));
        periods.put("Uut", new Integer(7));
        periods.put("Uuq", new Integer(7));
        periods.put("Uup", new Integer(7));
        periods.put("Uuh", new Integer(7));
        periods.put("Uus", new Integer(7));
        periods.put("Uuo", new Integer(7));
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
