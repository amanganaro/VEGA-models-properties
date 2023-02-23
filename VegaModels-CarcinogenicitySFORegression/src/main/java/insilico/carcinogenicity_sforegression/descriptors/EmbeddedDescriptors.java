package insilico.carcinogenicity_sforegression.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.carcinogenicity_sforegression.descriptors.weights.EState;
import insilico.carcinogenicity_sforegression.descriptors.weights.FunctionalGroups;
import insilico.carcinogenicity_sforegression.descriptors.weights.Mass;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.Arrays;


public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private final int MISSING_VALUE = -999;

    public double GATS2s = MISSING_VALUE;
    public double P_VSA_m_5 = MISSING_VALUE;
    public double ESpm8bo = MISSING_VALUE;
    public double nRNNOx = MISSING_VALUE;
    public double nFuranes = MISSING_VALUE;
    public double CATS2D_3_DA = MISSING_VALUE;
    public double CATS2D_4_AL = MISSING_VALUE;

    private final static String TYPE_A = "A";
    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";

    private double[][] ConnAugMatrix;


    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateFunctionalGroups(Mol);
        CalculateCats2D(Mol);
        CalculateGATS(Mol);
        CalculatePVSA(Mol);
        CalculateEAC(Mol);
    }

    private void CalculateCats2D(InsilicoMolecule Mol) {

        CATS2D_3_DA = 0; CATS2D_4_AL = 0;

        String[][] AtomCouplesAL= {
                {TYPE_A, TYPE_L},
        };

        String[][] AtomCouplesDA = {
                {TYPE_D, TYPE_A}
        };

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_3_DA = MISSING_VALUE; CATS2D_4_AL = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DA = MISSING_VALUE; CATS2D_4_AL = MISSING_VALUE;
            return;
        }
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DA = MISSING_VALUE; CATS2D_4_AL = MISSING_VALUE;
            return;
        }

        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (String[] strings : AtomCouplesDA) {

            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], strings[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (isIn(CatsTypes[j], strings[1])) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++)
                if (i == 3)
                    CATS2D_3_DA = desc[i];

        }

        for (String[] strings : AtomCouplesAL) {

            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], strings[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (isIn(CatsTypes[j], strings[1])) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++)
                if (i == 4)
                    CATS2D_4_AL = desc[i];

        }


    }

    private void CalculateEAC(InsilicoMolecule Mol){

        ESpm8bo = 0;
        final int MinSM = 1;
        final int MaxSM = 15;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ESpm8bo = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            ESpm8bo = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            ESpm8bo = MISSING_VALUE;
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeDegreeMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                if (EdgeAdjMat[i][j][0] != 0)
                    EdgeDegreeMat[i][j] = MoleculeUtilities.Bond2Double(m.getBond(j));

        DataMatrix = new Matrix(EdgeDegreeMat);

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        for (int i=MinSM; i<=MaxSM; i++) {
            double curSM = 0;
            for (int k=(eigenvalues.length-1); k>=0; k--) {
                curSM += Math.pow(eigenvalues[k], (i));
            }
            curSM = Math.log(1 + curSM);

            if ( i == 8)
                ESpm8bo = curSM;
        }
    }

    private void CalculatePVSA(InsilicoMolecule Mol){
        P_VSA_m_5 = 0;
        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            P_VSA_m_5 = MISSING_VALUE;
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

        double[] w = Mass.getWeights(m);

        int bins = BuildBinSize(0);

        for (int b=1; b<=bins; b++) {

            for (int i=0; i<nSK; i++) {
                if (w[i] == Descriptor.MISSING_VALUE) continue;
                if ( b == CalculateBin(0, w[i]))
                    if(b==5)
                        P_VSA_m_5 += VSA[i];
            }
        }
    }



    private void CalculateFunctionalGroups(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorBlock descriptorBlock = new FunctionalGroups();
        descriptorBlock.Calculate(Mol);
        nRNNOx = descriptorBlock.GetByName("nRNNOx").getValue();
        nFuranes = descriptorBlock.GetByName("nFuranes").getValue();
    }

    private void CalculateGATS(InsilicoMolecule Mol){

        GATS2s = 0;

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            GATS2s = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            GATS2s = MISSING_VALUE;
            return;
        }

        // !!! in origine era usata la topological matrix del cdk

        int nSK = m.getAtomCount();

        double[] w;
        try {
            EState ES = new EState(m);
            w = ES.getIS();

            // correction for compatibility with D7
            // H I-state is always 1
            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                    w[i] = 1;
            }
        } catch (Exception e) {
            w = new double[nSK];
            for (int i=0; i<nSK; i++) w[i]=Descriptor.MISSING_VALUE;
        }

        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;

        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 2) {
                    GearyAC += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                GearyAC = 0;
            } else {
                GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
            }
        }



        GATS2s = GearyAC;

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

    private int BuildBinSize(int curWeight) {
        if (curWeight == 0)
            return 5;
        if (curWeight == 1)
            return 8;
        if (curWeight == 2)
            return 4;
        if (curWeight == 3)
            return 8;
        return 0;
    }

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

                int nSglBnd = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrix[i][j]>0) {
                        if (ConnAugMatrix[i][j] == 1)
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
                    if (ConnAugMatrix[i][j]>0) {
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
                            if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
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
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] > 1.5) ) {
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
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] != 1) ) {
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
            if (tA) AtomTypes[i].add(TYPE_A);
            if (tD) AtomTypes[i].add(TYPE_D);
            if (tL) AtomTypes[i].add(TYPE_L);

        }

        return AtomTypes;
    }

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }

}
