package insilico.persistence_soil_quantitative_irfmn.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.persistence_soil_quantitative_irfmn.descriptors.weights.BurdenMatrix;
import insilico.persistence_soil_quantitative_irfmn.descriptors.weights.DescriptorMLogP;
import insilico.persistence_soil_quantitative_irfmn.descriptors.weights.IonizationPotential;
import insilico.persistence_soil_quantitative_irfmn.descriptors.weights.Polarizability;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double MLogP = MISSING_VALUE;
    private double SpPosA_p = MISSING_VALUE;
    private double B1_C_Cl = MISSING_VALUE;
    private double CATS2D_3_DL = MISSING_VALUE;
    private double P_VSA_i_3 = MISSING_VALUE;
    private double EEig7ri = MISSING_VALUE;
    private double B3_O_O = MISSING_VALUE;
    private double B1_C_O = MISSING_VALUE;

    private double[][] ConnAugMatrixCats2D;
    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";
    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateMLogP(Mol);
        CalculateTopoDistances(Mol);
        CalculateBE(Mol);
        CalculateCATS2D(Mol);
        CalculatePVSA(Mol);
        CalculateEAAC(Mol);
    }

    private void CalculateEAAC(InsilicoMolecule Mol){
        this.setEEig7ri(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setEEig7ri(MISSING_VALUE);
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            this.setEEig7ri(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setEEig7ri(MISSING_VALUE);
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeResMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeResMat.length; i++)
            for (int j=0; j<EdgeResMat[0].length; j++)
                EdgeResMat[i][j] = EdgeAdjMat[i][j][0]; // from standard edge adj matrix

        for (int i=0; i<m.getBondCount(); i++)
            EdgeResMat[i][i] = GetResonanceIntegral(m.getBond(i));

        DataMatrix = new Matrix(EdgeResMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        // EEig
        for (int i=1; i<=15; i++) {
            if(i == 7){
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    this.setEEig7ri(eigenvalues[idx]);
                else
                    this.setEEig7ri(0);
                return;
            }

        }
    }

    private void CalculatePVSA(InsilicoMolecule Mol){
        this.setP_VSA_i_3(0);

        int bins = 4;

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            this.setP_VSA_i_3(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();


        // Calculate VSA
        double[] VSA = new double[nSK];

        for (int i=0; i<nSK; i++) {

            IAtom at = m.getAtom(i);

            // R - van der waals radius
            double vdwR = GetVdWRadius(m, at);
            if (vdwR == Descriptor.MISSING_VALUE) continue;

            double coef = 0;
            for (IAtom connAt : m.getConnectedAtomsList(at)) {

                double connAt_vdwR = GetVdWRadius(m, connAt);
                if (connAt_vdwR == Descriptor.MISSING_VALUE) continue;

                // refR - reference bond length
                double refR = this.GetRefBondLength(at, connAt);
                if (refR == Descriptor.MISSING_VALUE) continue;

                // c - correction
                double c = 0;
                double bnd = MoleculeUtilities.Bond2Double(m.getBond(at, connAt));
                if (bnd == 1.5) c = 0.1;
                if (bnd == 2) c = 0.2;
                if (bnd == 3) c = 0.3;

                double g_1 = Math.max(  Math.abs(vdwR - connAt_vdwR) , GetRefBondLength(at, connAt) - c ) ;
                double g_2 = vdwR + connAt_vdwR;
                double g  = Math.min(g_1, g_2) ;

                coef += ( Math.pow(connAt_vdwR,2) - Math.pow( (vdwR - g), 2) ) / (g) ;
            }

            VSA[i] = ( 4.0 * Math.PI * Math.pow(vdwR,2) ) - Math.PI * vdwR * coef ;
        }

        double[] w = IonizationPotential.getWeightsNormalized(m);

        for (int b=1; b<=bins; b++) {

            if (b==3) {
                double PVSA = 0;
                for (int i=0; i<nSK; i++) {
                    if (w[i] == MISSING_VALUE) continue;
                    if ( b == CalculateBin(2, w[i]))
                        PVSA += VSA[i];
                }
                this.setP_VSA_i_3(PVSA);
                return;
            }
        }

    }

    private void CalculateCATS2D(InsilicoMolecule Mol){
        this.setCATS2D_3_DL(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }
        try {
            ConnAugMatrixCats2D = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (int d=0; d<AtomCouples.length; d++) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i=0; i<nSK; i++) {
                if ( isIn(CatsTypes[i], AtomCouples[d][0]) ) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if ( isIn(CatsTypes[j], AtomCouples[d][1]) ) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i=0; i<desc.length; i++)
                if(i == 3)
                    this.setCATS2D_3_DL(desc[i]);

        }
    }

    private void CalculateBE(InsilicoMolecule Mol){
        this.setSpPosA_p(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setSpPosA_p(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            this.setSpPosA_p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = Polarizability.getWeights(m);

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        double posSum = 0;
        for (double e : eigenvalues) {
            if (e>0) posSum += e;
        }

        double SpPosA = posSum / nSK;
        this.setSpPosA_p(SpPosA);


    }

    private void CalculateTopoDistances(InsilicoMolecule Mol){
        this.setB1_C_Cl(0);
        this.setB3_O_O(0);
        this.setB1_C_O(0);

        String[][] AtomCouplesCO = {
                {"C", "O"},
        };

        String[][] AtomCouplesOO = {
                {"O", "O"},
        };

        String[][] AtomCouplesCCl= {
                {"C", "Cl"},
        };


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setB1_C_Cl(MISSING_VALUE);
            this.setB3_O_O(MISSING_VALUE);
            this.setB1_C_O(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setB1_C_Cl(MISSING_VALUE);
            this.setB3_O_O(MISSING_VALUE);
            this.setB1_C_O(MISSING_VALUE);
            return;
        }

        int[] descB = new int[10];
        Arrays.fill(descB, 0);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesCO[0][0])) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesCO[0][1])) {

                        // B (presence of pair) and F (number of couples)
                        if (TopoMat[i][j] <= 10) {
                            descB[TopoMat[i][j]-1] = 1;
                        }

                    }
                }
            }
        }


        for (int i=0; i<descB.length; i++) {
            if(i == 0)
                this.setB1_C_O(descB[i]);
        }

        descB = new int[10];
        Arrays.fill(descB, 0);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesOO[0][0])) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesOO[0][1])) {

                        // B (presence of pair) and F (number of couples)
                        if (TopoMat[i][j] <= 10) {
                            descB[TopoMat[i][j]-1] = 1;
                        }

                    }
                }
            }
        }


        for (int i=0; i<descB.length; i++) {
            if(i == 2)
                this.setB3_O_O(descB[i]);
        }

        descB = new int[10];
        Arrays.fill(descB, 0);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesCCl[0][0])) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesCCl[0][1])) {

                        // B (presence of pair) and F (number of couples)
                        if (TopoMat[i][j] <= 10) {
                            descB[TopoMat[i][j]-1] = 1;
                        }

                    }
                }
            }
        }


        for (int i=0; i<descB.length; i++) {
            if(i == 0)
                this.setB1_C_Cl(descB[i]);
        }

    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        try {
            this.MLogP = 0;
            DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
            this.setMLogP(descriptorMLogP.getMLogP());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setMLogP(MISSING_VALUE);
        }

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
    private ArrayList<String>[] setCatsAtomType(IAtomContainer m) {

        int nSK = m.getAtomCount();
        ArrayList<String>[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false;

            // Definition of CATS types
            //
            // A: O, N without H
            // N: [+], NH2
            // P: [-], COOH, POOH, SOOH

            // Hydrogens
            int H = 0;
            try {
                H = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) { }

            // [+]
            if (CurAt.getFormalCharge() > 0) {

                boolean NpOm = false;
                if (ConnAugMatrixCats2D[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrixCats2D[i][j]==1) {
                            if (ConnAugMatrixCats2D[j][j] == 8) {
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

                int nSglBnd = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if (ConnAugMatrixCats2D[i][j] == 1)
                            nSglBnd++;
                        else
                            nOtherBnd++;
                    }
                }

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (CurAt.getFormalCharge() == 0) &&( (H == 1) || (H ==2) ) )
                    tD = true;

            }

            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                    (CurAt.getFormalCharge() == 0) )  {

                int nSglBnd = 0, nDblO = 0, nSglOH = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if (ConnAugMatrixCats2D[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrixCats2D[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrixCats2D[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrixCats2D[i][j] == 2) && (ConnAugMatrixCats2D[j][j] == 8) )
                                nDblO++;
                            else
                                nOtherBnd++;
                        }
                    }
                }

                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherBnd == 0) )
                    tN = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                boolean connOnlyToSingleC = true;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if ( (ConnAugMatrixCats2D[j][j] != 6) || (ConnAugMatrixCats2D[i][j] > 1.5) ) {
                            connOnlyToSingleC = false;
                            break;
                        }
                    }
                }
                if (connOnlyToSingleC)
                    tL = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                boolean connOnlyToSingleC = true;
                int nSingleC = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if ( (ConnAugMatrixCats2D[j][j] != 6) || (ConnAugMatrixCats2D[i][j] != 1) ) {
                            connOnlyToSingleC = false;
                            break;
                        } else {
                            nSingleC++;
                        }
                    }
                }
                if ( (connOnlyToSingleC) && (nSingleC == 2) )
                    tL = true;
            }


            // Sets final types
            if (tD) AtomTypes[i].add(TYPE_D);
            if (tL) AtomTypes[i].add(TYPE_L);

        }

        return AtomTypes;
    }

    private int CalculateBin(int curWeight, double value) {

        if (curWeight == 0) {
            // Mass
            if (value < 1) return 1;
            if (value < 1.2) return 2;
            if (value < 1.6) return 3;
            if (value < 3) return 4;
            return 5;
        }

        if (curWeight == 1) {
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

        if (curWeight == 2) {
            // Ionization potential
            if (value < 1) return 1;
            if (value < 1.15) return 2;
            if (value < 1.25) return 3;
            return 4;
        }

        if (curWeight == 3) {
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

        return 0;
    }


    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }

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
            if (Manipulator.CountImplicitHydrogens(at) > 0)
                return 2.152;

            return 1.779;
        }

        if (s.equalsIgnoreCase("H")) {
            IAtom connAt = m.getConnectedAtomsList(at).get(0);
            if (connAt.getSymbol().equalsIgnoreCase("0"))
                return 0.8;
            if (connAt.getSymbol().equalsIgnoreCase("N"))
                return 0.7;
            if (connAt.getSymbol().equalsIgnoreCase("P"))
                return 0.7;
            return 1.485;
        }

        return Descriptor.MISSING_VALUE;
    }

    private double GetResonanceIntegral(IBond bnd) {

        IAtom atA =  bnd.getAtom(0);
        IAtom atB =  bnd.getAtom(1);
        String A = atA.getSymbol();
        String B = atB.getSymbol();

        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("C") == 0)) )
            return 1.00;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("B") == 0)) ||
                ((A.compareToIgnoreCase("B") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("N") == 0)) ||
                ((A.compareToIgnoreCase("N") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.9;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("O") == 0)) ||
                ((A.compareToIgnoreCase("O") == 0) && (B.compareToIgnoreCase("C") == 0))) {
            if (bnd.getOrder() == IBond.Order.SINGLE)
                return 0.8;
            if (bnd.getOrder() == IBond.Order.DOUBLE)
                return 1.2;
        }
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("S") == 0)) ||
                ((A.compareToIgnoreCase("S") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("F") == 0)) ||
                ((A.compareToIgnoreCase("F") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Cl") == 0)) ||
                ((A.compareToIgnoreCase("Cl") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.4;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Br") == 0)) ||
                ((A.compareToIgnoreCase("Br") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.3;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("I") == 0)) ||
                ((A.compareToIgnoreCase("I") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.1;

        return 0.00;
    }




}
