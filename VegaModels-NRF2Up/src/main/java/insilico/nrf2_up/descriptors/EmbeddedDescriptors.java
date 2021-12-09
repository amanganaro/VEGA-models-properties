package insilico.nrf2_up.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.descriptor.blocks.*;
import insilico.descriptor.blocks.logP.MLogP;
import insilico.descriptor.localization.StringSelectorDescriptors;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private final double MISSING_VALUE = -999;

    public double MLOGP = MISSING_VALUE;
    public double P_VSA_ppp_L = MISSING_VALUE;
    public double P_VSA_e_2 = MISSING_VALUE;
    public double SpMin2_Bh_e= MISSING_VALUE;
    public double P_VSA_i_2 = MISSING_VALUE;
    public double F07_C_C = MISSING_VALUE;
    public double NaasC = MISSING_VALUE;
    public double IDDE = MISSING_VALUE;
    public double CATS2D_06_LL = MISSING_VALUE;
    public double CATS2D_04_LL = MISSING_VALUE;
    public double IC3 = MISSING_VALUE;
    public double ATSC2s = MISSING_VALUE;
    public double P_VSA_ppp_cyc = MISSING_VALUE;

    public double[] getDescriptors(){
        return new double[]{MLOGP, P_VSA_ppp_L, P_VSA_e_2, SpMin2_Bh_e, P_VSA_i_2, F07_C_C, NaasC, IDDE, CATS2D_06_LL, CATS2D_04_LL, IC3, ATSC2s, P_VSA_ppp_cyc};
    }

    public EmbeddedDescriptors(InsilicoMolecule mol,boolean fromFile) throws MalformedURLException {
        if(fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
//        System.out.println();
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateMLogP(mol);
        CalculatePVSA(mol);
        CalculateSpMin(mol);
        CalculateF(mol);
        CalculateNaasc(mol);
        CalculateIC(mol);
        CalculateATSC(mol);
        CalculateCATS(mol);
    }

    private void CalculateCATS(InsilicoMolecule mol) {

        String[] TYPE_L = { "L", ""};
        String[][][] AtomCouples = {{TYPE_L, TYPE_L}};
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        int nSK = curMol.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        double[][] ConnAugMatrix = null;
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        ArrayList<String>[] CatsTypes = setCatsAtomType(curMol, ConnAugMatrix);


        for (String[][] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], atomCouple[0][0])) {
                    for (int j = i; j < nSK; j++) {
//                        if (i==j) continue;
                        if (isIn(CatsTypes[j], atomCouple[1][0])) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++){
                if (i == 4)
                    CATS2D_04_LL = desc[i];
                if(i==6)
                    CATS2D_06_LL = desc[i];
            }
        }

    }

    private void CalculateMLogP(InsilicoMolecule mol) {
        DescriptorBlock block = new MLogP();
        try {
            block.Calculate(mol);
            MLOGP = block.GetByName("MLogP").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculatePVSA(InsilicoMolecule mol) {

        // P_VSA are calculated on H filled molecules
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        int nSK = m.getAtomCount();
        int nBO = m.getBondCount();

        GhoseCrippenACF GC = new GhoseCrippenACF(m, true);
        int[] ACF = GC.GetACF();

        double[] VSA = new double[nSK];
        double[] PartialVSA = new double[nSK];
        double[] IdealBondLen = new double[nBO];
        double[] VdWRadius = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = m.getAtom(i);

            VSA[i] = 0;
            PartialVSA[i] = 0;

            // R - van der waals radius
            VdWRadius[i] = GetVdWRadius(m, at);
            if (VdWRadius[i] == Descriptor.MISSING_VALUE) {
                log.warn("Missing VSA");
                VSA[i] = Descriptor.MISSING_VALUE;
                PartialVSA[i] = Descriptor.MISSING_VALUE;
            }
        }

        for (int i=0; i<nBO; i++) {
            IBond bo = m.getBond(i);

            // refR - reference bond length
            double refR = this.GetRefBondLength(bo.getAtom(0), bo.getAtom(1));
            if (refR == Descriptor.MISSING_VALUE) {
                log.warn("Missing radius");
            }

            // c - correction
            double c = 0;
            if (refR != 0) {
                double bnd = MoleculeUtilities.Bond2Double(bo);
                if (bnd == 1.5) c = 0.1;
                if (bnd == 2) c = 0.2;
                if (bnd == 3) c = 0.3;
            }

            IdealBondLen[i] = refR - c;
        }

        for (int b=0; b<nBO; b++) {
            IBond bo = m.getBond(b);


            int at1 = m.indexOf(bo.getAtom(0));
            int at2 = m.indexOf(bo.getAtom(1));
            if ( (IdealBondLen[b] == 0) || (VdWRadius[at1] == Descriptor.MISSING_VALUE) || (VdWRadius[at2] == Descriptor.MISSING_VALUE) ) {
                PartialVSA[at1] = Descriptor.MISSING_VALUE;
                PartialVSA[at2] = Descriptor.MISSING_VALUE;
                VSA[at1] = Descriptor.MISSING_VALUE;
                VSA[at2] = Descriptor.MISSING_VALUE;
                continue;
            }

            double curDistance;
            double diffRadius = Math.abs(VdWRadius[at1]-VdWRadius[at2]);
            if (diffRadius > IdealBondLen[b])
                curDistance = diffRadius;
            else
                curDistance = IdealBondLen[b];

            if (curDistance > (VdWRadius[at1]+VdWRadius[at2]))
                curDistance = VdWRadius[at1]+VdWRadius[at2];

            PartialVSA[at1] = PartialVSA[at1] + ( Math.pow(VdWRadius[at2],2) - Math.pow((VdWRadius[at1] - curDistance), 2) ) / curDistance;
            PartialVSA[at2] = PartialVSA[at2] + ( Math.pow(VdWRadius[at1],2) - Math.pow((VdWRadius[at2] - curDistance), 2) ) / curDistance;
        }

        for (int i=0; i<nSK; i++) {
            if (VSA[i] != Descriptor.MISSING_VALUE)
                VSA[i] = 4 * Math.PI * Math.pow(VdWRadius[i], 2) - (Math.PI * VdWRadius[i] * PartialVSA[i]);
            if (VSA[i] < 0)
                VSA[i] = 0;
        }

        // weight i
//        P_VSA_i_2
        double[] w = new WeightsIonizationPotential().getScaledWeights(m);
        int bins = 4;
        for (int b=1; b<=bins; b++) {
            double PVSA = 0;
            for (int i = 0; i < nSK; i++) {
                if (w[i] == Descriptor.MISSING_VALUE) continue;
                if (b == CalculateBin("i", w[i])) {
                    PVSA += VSA[i];
                    if(b == 2)
                        P_VSA_i_2 = PVSA;
                }
            }
        }

        //
//        P_VSA_e_2
        w = new WeightsElectronegativity().getScaledWeights(m);
        bins = 6;
        for (int b=1; b<=bins; b++) {
            double PVSA = 0;
            for (int i = 0; i < nSK; i++) {
                if (w[i] == Descriptor.MISSING_VALUE) continue;
                if (b == CalculateBin("i", w[i])) {
                    PVSA += VSA[i];
                    if(b == 2)
                        P_VSA_e_2 = PVSA;
                }
            }
        }

        IAtomContainer curMolNoH;
        double[][] ConnAugMatrixNoH = null;
        try {
            curMolNoH = mol.GetStructure();
            ConnAugMatrixNoH = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        int nSKnoH = curMolNoH.getAtomCount();

        Cats2D cats = new Cats2D();
        ArrayList<String>[] AtomTypesOnMolWithoutH = cats.setCatsAtomType(curMolNoH, ConnAugMatrixNoH);

        String curAtomType = "L";
        double PVSA = 0;

        for (int i = 0; i < nSKnoH; i++) {

            for (String atType : AtomTypesOnMolWithoutH[i]){
                if (atType.equalsIgnoreCase(curAtomType)) {
                    PVSA += VSA[i];
                    break;
                }
            }
        }

        P_VSA_ppp_L = PVSA;

        curAtomType = "cyc";
        PVSA = 0;

        for (int i = 0; i < nSKnoH; i++) {

            for (String atType : AtomTypesOnMolWithoutH[i]){
                if (atType.equalsIgnoreCase(curAtomType)) {
                    PVSA += VSA[i];
                    break;
                }
            }
        }

        P_VSA_ppp_cyc = PVSA;
    }

    private void CalculateSpMin(InsilicoMolecule mol) {
//        DescriptorBlock block = new BurdenEigenvalues();
//        try {
//            block.Calculate(mol);
//            SpMin2_Bh_e = block.GetByName("SpMin2_Bh(e)").getValue();
//        } catch (DescriptorNotFoundException ex){
//            log.warn(ex.getMessage());
//        }
        IAtomContainer m = null;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
        }

        // Gets matrix
        double[][] BurdenMat = new double[0][];
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
        }

        int nSK = m.getAtomCount();
        iBasicWeight weight = new WeightsElectronegativity();
        double[] w = weight.getScaledWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH, valL;
            if (i > (eigenvalues.length - 1)) {
                valH = 0;
                valL = 0;
            } else {
                if (eigenvalues[eigenvalues.length - 1 - i] > 0)
                    valH = eigenvalues[eigenvalues.length - 1 - i];
                else
                    valH = 0;
                if (eigenvalues[i] < 0)
                    valL = Math.abs(eigenvalues[i]);
                else
                    valL = 0;
            }

            if (i == 1)
                SpMin2_Bh_e = valL;

        }

        //TODO viene diverso da Dragon!

    }

    private void CalculateF(InsilicoMolecule mol) {
        DescriptorBlock block = new AtomPairs2D();
        try {
            block.Calculate(mol);
            F07_C_C = block.GetByName("F07[C-C]").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        String[][] ATOM_COUPLES = {
                {"C", "C"}
        };

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn("Invalid structure, unable to calculate: " + mol.GetSMILES());
        }
        int MAX_TOPO_DISTANCE = 10;
        for (int d = 0; d< ATOM_COUPLES.length; d++) {

            int descT = 0;
            int[] descB = new int[MAX_TOPO_DISTANCE];
            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2) // DA VEDERE PERCHE MAGGIORE DI 2
                                descT += TopoMat[i][j];

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j] - 1] = 1;
                                descF[TopoMat[i][j] - 1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (ATOM_COUPLES[d][0].compareTo(ATOM_COUPLES[d][1]) == 0) {
                descT /= 2;
                for (int i = 0; i < descF.length; i++)
                    descF[i] /= 2;
            }
            for (int i=0; i<descB.length; i++) {
                int lag = i + 1;
                if (lag == 7)
                    F07_C_C = descF[i];
            }
        }
    }

    private void CalculateNaasc(InsilicoMolecule mol) {
        DescriptorBlock block = new EStateIndices();
        try {
            block.Calculate(mol);
            NaasC = block.GetByName("NaasC").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateIC(InsilicoMolecule mol) {
        DescriptorBlock block = new InformationContent();
        try {
            block.Calculate(mol);
            IDDE = block.GetByName("IDDE").getValue();
            IC3 = block.GetByName("IC3").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        // Gets matrices
        int[][] TopoDistMat;
        try {
            TopoDistMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = curMol.getAtomCount();
        int[] VertexDeg = WeightsVertexDegree.getWeights(curMol, false);
        double[] IStates = (new WeightsIState()).getWeights(curMol, false);

        // rGes
        double[] ETIS = new double[nSK];
        for (int i=0; i<nSK; i++) {
            ETIS[i] = IStates[i];
            if (ETIS[i] != Descriptor.MISSING_VALUE) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    ETIS[i] += (IStates[i] - IStates[j]) / Math.pow(TopoDistMat[i][j] + 1, 2);
                }
            }
        }

        double[] s0k_Equal = new double[nSK];
        int[] ETISequal = new int[nSK];
        for (int i=0; i<nSK; i++)
            s0k_Equal[i] = ETIS[i];
        int s0k_Cequiv = 0;

        for (int i=0; i<nSK; i++) {
            if (s0k_Equal[i] != Descriptor.MISSING_VALUE) {
                int k = 1;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (Math.abs(s0k_Equal[i] - s0k_Equal[j]) < 0.001) {
                        k++;
                        s0k_Equal[j] = Descriptor.MISSING_VALUE;
                    }
                }
                ETISequal[s0k_Cequiv] = k;
                s0k_Equal[i] = Descriptor.MISSING_VALUE;
                s0k_Cequiv++;
            }
        }


        // IDDE
        int[] DistDeg = new int[nSK];
        for (int i=0; i<nSK; i++) {
            DistDeg[i] = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if (TopoDistMat[i][j] > 0)
                    DistDeg[i] += TopoDistMat[i][j];
            }
        }

        int[] Equal = new int[nSK];
        int[] DistDegEq = new int[nSK];
        for (int i=0; i<nSK; i++)
            Equal[i] = DistDeg[i];

        for (int i=0; i<nSK; i++) {
            if (Equal[i] != 0) {
                int k = 1;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (Equal[i] == Equal[j]) {
                        k++;
                        Equal[j] = 0;
                    }
                }
                DistDegEq[i] = k;
                Equal[i] = 0;
            }
        }

        double IDDE = 0;
        for (int i=0; i<nSK; i++) {
            if (DistDegEq[i] > 0)
                IDDE -= ((double) DistDegEq[i] / (double) nSK) * Math.log(((double) DistDegEq[i] / (double) nSK));
        }
        this.IDDE = 1.0 / Math.log(2) * IDDE;


        IAtomContainer curMolH;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            curMolH = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        // Gets matrices
        double[][] ConnMatH;
        int[][] TopoDistMatH;
        try {
            ConnMatH= ConnectionAugMatrix.getMatrix(curMolH);
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSKH = curMolH.getAtomCount();
        int[] VertexDegH = WeightsVertexDegree.getWeights(curMolH, false);

        double bic_denom = 0;
        for (IBond bnd : curMolH.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag=1; CurLag<=5; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSKH);
            for (int i=0; i<nSKH; i++) {
                IAtom atStart =  curMolH.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j=0; j<nSKH; j++) {
                    if (i==j) continue;
                    if (TopoDistMatH[i][j] == CurLag) {
                        IAtom atEnd =  curMolH.getAtom(j);
                        ShortestPaths shortestPaths = new ShortestPaths(curMolH, atStart);
                        List<IAtom> sp = Arrays.asList(shortestPaths.atomsTo(atEnd));
                        StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                        for (int k=0; k<(sp.size()-1); k++) {
                            int a = curMolH.indexOf(sp.get(k));
                            int b = curMolH.indexOf(sp.get(k + 1));
                            if (ConnMatH[a][b] == 1)
                                bufPath.append("s");
                            if (ConnMatH[a][b] == 2)
                                bufPath.append("d");
                            if (ConnMatH[a][b] == 3)
                                bufPath.append("t");
                            if (ConnMatH[a][b] == 1.5)
                                bufPath.append("a");
                            bufPath.append(sp.get(k + 1).getAtomicNumber());
                            bufPath.append("(").append(VertexDegH[curMolH.indexOf(sp.get(k + 1))]).append(")");
                        }
                        CurNeig.add(bufPath.toString());
                    }
                }
                Collections.sort(CurNeig);
                NeigList.add(CurNeig);
            }

            // Calculates equivalence classes
            ArrayList<ArrayList<String>> G = new ArrayList<>();
            ArrayList<Integer> Gn = new ArrayList<>();
            for (int i=0; i<nSKH; i++) {
                ArrayList<String> CurNeig = NeigList.get(i);
                boolean foundMatch = false;
                for (int k=0; k<G.size(); k++) {
                    if (CompareNeigVector(CurNeig, G.get(k))) {
                        foundMatch = true;
                        int buf = Gn.get(k);
                        Gn.set(k, (buf+1));
                        break;
                    }
                }
                if (!foundMatch) {
                    G.add(CurNeig);
                    Gn.add(1);
                }
            }

            // Calculate IC and CIC indices

            double ic=0;
            for (int i=0; i<Gn.size(); i++)
                ic += ((double)Gn.get(i)/nSKH) * (Math.log((double)Gn.get(i)/nSKH));
            ic = (-1.00 / Math.log(2)) * ic;

            if(CurLag == 3)
                IC3 = ic;
        }

    }

    private void CalculateATSC(InsilicoMolecule mol) {
//        DescriptorBlock block = new AutoCorrelation();
//        try {
//            block.Calculate(mol);
//            ATSC2s = block.GetByName("ATSC2s").getValue();
//        } catch (DescriptorNotFoundException ex){
//            log.warn(ex.getMessage());
//        }

        try {



            // AC are calculated on H-filled molecules

            IAtomContainer m;
            try {
                IAtomContainer orig_m = mol.GetStructure();
                m = Manipulator.AddHydrogens(orig_m);
            } catch (InvalidMoleculeException | GenericFailureException e) {
                log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
                return;
            }

            // Gets matrices
            int[][] TopoMatrix;
            try {
                TopoMatrix = TopoDistanceMatrix.getMatrix(m);
            } catch (Exception e) {
                log.warn(e.getMessage());
                return;
            }

            int nSK = m.getAtomCount();

            // Cycle for all found weighting schemes


            // Sets needed weights
            double[] w = new WeightsIState().getWeights(m, true);



            // correction for compatibility with D7
            // H I-state is always 1
            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                    w[i] = 1;
            }

            // Calculates weights averages
            double wA = 0;
            for (int i=0; i<nSK; i++)
                wA += w[i];
            wA = wA / ((double) nSK);

            // Calculates autocorrelations
//                for (int lag=1; lag<=8; lag++) {

            int lag = 2;

            double AC=0, ACS=0, MoranAC=0, GearyAC=0;
            double denom = 0, delta = 0;

            for (int i=0; i<nSK; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        AC += w[i] * w[j];
                        ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                        MoranAC += (w[i] - wA) * (w[j] - wA);
                        GearyAC += Math.pow((w[i] - w[j]), 2);
                        delta++;
                    }
            }

            if (delta > 0) {
                if (denom == 0) {
                    MoranAC = 1;
                    GearyAC = 0;
                } else {
                    MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                    GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
                }
            }

            // AC transformed in log form
            AC /= 2.0;
            AC = Math.log(1 + AC);

            ACS /= 2.0;

            ATSC2s = ACS;

        } catch (Throwable e) {
            log.warn("Unable to calculate: ATSC - " + e.getMessage());
        }


    }


    private void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {
        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-NRF2Up/src/main/resources/data/dataset.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            br.readLine();
            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[1].trim()).GetSMILES())) {

                    MLOGP = Double.parseDouble(lineArray[2]);
                    P_VSA_ppp_L = Double.parseDouble(lineArray[3]);
                    P_VSA_e_2 = Double.parseDouble(lineArray[4]);
                    SpMin2_Bh_e = Double.parseDouble(lineArray[5]);
                    P_VSA_i_2 = Double.parseDouble(lineArray[6]);
                    F07_C_C = Double.parseDouble(lineArray[7]);
                    NaasC = Double.parseDouble(lineArray[8]);
                    IDDE = Double.parseDouble(lineArray[9]);
                    CATS2D_06_LL = Double.parseDouble(lineArray[10]);
                    CATS2D_04_LL = Double.parseDouble(lineArray[11]);
                    IC3 = Double.parseDouble(lineArray[12]);
                    ATSC2s = Double.parseDouble(lineArray[13]);
                    P_VSA_ppp_cyc = Double.parseDouble(lineArray[14]);
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }


    }

    private boolean CompareNeigVector(ArrayList<String> A, ArrayList<String> B) {

        // note: integer in the vectors should be already sorted

        if (A.size()!=B.size())
            return false;
        else {
            for (int i=0; i<A.size(); i++)
                if (!A.get(i).equalsIgnoreCase(B.get(i)))
                    return false;
        }

        return true;
    }

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }


    /**
     * Sets CATS 2D atom types for each atom, as a list of string containing
     * all matching types for each atom.
     **/
    public ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        ArrayList[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false, tCyc=false;

            // Definition of CATS types
            //
            // A: O, N without H
            // N: [+], NH2
            // P: [-], COOH, POOH, SOOH

            // Hydrogens
            int H = 0;
            try {
                H = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn(StringSelectorDescriptors.getString("unable_count_h"));
            }

            // counters
            int nSglBnd = 0, nOtherBnd = 0, VD = 0;
            int nC = 0, nDblO = 0, nOtherNonOBond=0, nSglOH = 0;
            for (int j=0; j<nSK; j++) {
                if (j==i) continue;
                if (ConnAugMatrix[i][j]>0) {

                    VD++;

                    if (ConnAugMatrix[j][j] == 6)
                        nC++;

                    if (ConnAugMatrix[i][j] == 1) {
                        nSglBnd++;

                        if (ConnAugMatrix[j][j] == 8) {
                            int Obonds = 0;
                            for (int k=0; k<nSK; k++) {
                                if (k == j) continue;
                                if (ConnAugMatrix[k][j]>0) Obonds++;
                            }
                            if (Obonds == 1) nSglOH++;
                        }

                    } else {
                        nOtherBnd++;
                        if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
                            nDblO++;
                        else
                            nOtherNonOBond++;
                    }
                }

            }


            // [+]
            if (CurAt.getFormalCharge() > 0) {

                boolean NpOm = false;
                if (ConnAugMatrix[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrix[i][j]==1) {
                            if (ConnAugMatrix[j][j] == 8) {
                                IAtom Oxy = m.getAtom(j);
                                if (Oxy.getFormalCharge()!=0)
                                    NpOm = true;
                            }
                        }
                    }
                }

                if (!NpOm)
                    tP = true;
            }

            // [-]
            if (CurAt.getFormalCharge() < 0)
                tN = true;

            // O
            if (CurAt.getSymbol().equalsIgnoreCase("O")) {
                tA = true;

                if ( (CurAt.getFormalCharge() == 0) && (H == 1))
                    tD = true;

            }

            // N (NH2 and N without H)
            if (CurAt.getSymbol().equalsIgnoreCase("N")) {

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (H == 1) || (H ==2) )
                    tD = true;

            }

            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                    (CurAt.getFormalCharge() == 0) )  {

                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherNonOBond == 0) )
                    tN = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                if ( VD == nC)
                    tL = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                if ( (nC == 2) && ( (VD+H) == 2 ))
                    tL = true;
            }

            if (CurAt.isInRing())
                tCyc = true;


            // Sets final types
            if (tL) AtomTypes[i].add("L");

        }

        return AtomTypes;
    }

    private final static String[][] PPP_TYPES = {
            Cats2D.TYPE_D,
            Cats2D.TYPE_A,
            Cats2D.TYPE_P,
            Cats2D.TYPE_N,
            Cats2D.TYPE_L,
            Cats2D.TYPE_CYC
    };

    private final static Object[][] RefBondLengths = {
            {"Br","Br",2.54},
            {"Br","C",1.97},
            {"Br","Cl",2.36},
            {"Br","F",1.85},
            {"Br","H",1.44},
            {"Br","I",2.65},
            {"Br","N",1.84},
            {"Br","O",1.58},
            {"Br","P",2.37},
            {"Br","S",2.21},
            {"C","C",1.54},
            {"C","Cl",1.8},
            {"C","F",1.35},
            {"C","H",1.06},
            {"C","I",2.12},
            {"C","N",1.47},
            {"C","O",1.43},
            {"C","P",1.85},
            {"C","S",1.81},
            {"Cl","Cl",2.31},
            {"Cl","F",1.63},
            {"Cl","H",1.22},
            {"Cl","I",2.56},
            {"Cl","N",1.74},
            {"Cl","O",1.41},
            {"Cl","P",2.01},
            {"Cl","S",2.07},
            {"F","F",1.28},
            {"F","H",0.87},
            {"F","I",2.04},
            {"F","N",1.41},
            {"F","O",1.32},
            {"F","P",1.5},
            {"F","S",1.64},
            {"H","I",1.63},
            {"H","N",1.01},
            {"H","O",0.97},
            {"H","P",1.41},
            {"H","S",1.31},
            {"I","I",2.92},
            {"I","N",2.26},
            {"I","O",2.14},
            {"I","P",2.49},
            {"I","S",2.69},
            {"N","N",1.45},
            {"N","O",1.46},
            {"N","P",1.6},
            {"N","S",1.76},
            {"O","O",1.47},
            {"O","P",1.57},
            {"O","S",1.57},
            {"P","P",2.26},
            {"P","S",2.07},
            {"S","S",2.05}
    };

    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }


    private double GetRefBondLength(IAtom at1, IAtom at2) {
        double len = Descriptor.MISSING_VALUE;
        for (int i=0; i<RefBondLengths.length; i++) {
            if (AtomCouple(at1, at2, (String)RefBondLengths[i][0], (String)RefBondLengths[i][1])) {
                len = (Double)RefBondLengths[i][2];
                break;
            }
        }
        return len;
    }


    private double GetVdWRadius(IAtomContainer m, IAtom at) {

        String s = at.getSymbol();
        if (s.equalsIgnoreCase("C"))
            return 1.950;
        if (s.equalsIgnoreCase("N"))
            return 1.950;
        if (s.equalsIgnoreCase("F"))
            return  1.496;
        if (s.equalsIgnoreCase("P"))
            return  2.287;
        if (s.equalsIgnoreCase("S"))
            return  2.185;
        if (s.equalsIgnoreCase("Cl"))
            return  2.044;
        if (s.equalsIgnoreCase("Br"))
            return  2.166;
        if (s.equalsIgnoreCase("I"))
            return 2.358;

        if (s.equalsIgnoreCase("O")) {

            // oxide
            if (m.getConnectedAtomsList(at).size() == 1)
                return 1.810;

            // acid
            int nH = 0;
            for (IAtom con : m.getConnectedAtomsList(at)) {
                if (con.getSymbol().equalsIgnoreCase("H"))
                    nH++;
            }
            if (nH > 0)
                return 2.152;

            return 1.779;
        }

        if (s.equalsIgnoreCase("H")) {
            IAtom connAt = m.getConnectedAtomsList(at).get(0);
            if (connAt.getSymbol().equalsIgnoreCase("O"))
                return 0.8;
            if (connAt.getSymbol().equalsIgnoreCase("N"))
                return 0.7;
            if (connAt.getSymbol().equalsIgnoreCase("P"))
                return 0.7;
            return 1.485;
        }

        return Descriptor.MISSING_VALUE;
    }

    private int CalculateBin(String curWeight, double value) {

        if (curWeight.equalsIgnoreCase("LogP")) {
            // LOGP
            if (value < -1.5) return 1;
            if (value < -0.5) return 2;
            if (value < -0.25) return 3;
            if (value < 0) return 4;
            if (value < 0.25) return 5;
            if (value < 0.52) return 6;
            if (value < 0.75) return 7;
            return 8;
        }

        if (curWeight.equalsIgnoreCase("MR")) {
            // Molar Refractivity
            if (value < 0.9) return 1;
            if (value < 1.5) return 2;
            if (value < 2.0) return 3;
            if (value < 2.5) return 4;
            if (value < 3.0) return 5;
            if (value < 4.0) return 6;
            if (value < 6.0) return 7;
            return 8;
        }

        if (curWeight.equalsIgnoreCase("m")) {
            // Mass
            if (value < 1.0) return 1;
            if (value < 1.2) return 2;
            if (value < 1.6) return 3;
            if (value < 3.0) return 4;
            return 5;
        }

        if (curWeight.equalsIgnoreCase("v")) {
            // VdW volume
            if (value < 0.5) return 1;
            if (value < 1.0) return 2;
            if (value < 1.3) return 3;
            return 4;
        }

        if (curWeight.equalsIgnoreCase("e")) {
            // Sanderson electronegativity
            if (value < 1.0) return 1;
            if (value < 1.1) return 2;
            if (value < 1.2) return 3;
            if (value < 1.3) return 4;
            if (value < 1.4) return 5;
            return 6;
        }

        if (curWeight.equalsIgnoreCase("p")) {
            // Polarizability
            if (value < 0.4) return 1;
            if (value < 1.0) return 2;
            if (value < 2.0) return 3;
            return 4;
        }

        if (curWeight.equalsIgnoreCase("i")) {
            // Ionization potential
            if (value < 1) return 1;
            if (value < 1.15) return 2;
            if (value < 1.25) return 3;
            return 4;
        }

        if (curWeight.equalsIgnoreCase("s")) {
            // I-State
            if (value < 1.0) return 1;
            if (value < 1.5) return 2;
            if (value < 2.0) return 3;
            if (value < 2.5) return 4;
            if (value < 3.0) return 5;
            return 6;
        }

        return 0;
    }



}
