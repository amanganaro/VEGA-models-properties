package insilico.algae_ec50.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.algae_ec50.descriptors.weights.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    private double[][] ConnAugMatrix;

    public double ALogP;
    public double MLogP;
    public double BEH5p;
    public double MATS4s;
    public double BEL5e;
    public double MATS3s;
    public double X3v;
    public double Mw;
    public double CATS2D_3_DL;
    public double P_VSA_m_4;

    public EmbeddedDescriptors() {
        ALogP = MISSING_VALUE;
        MLogP = MISSING_VALUE;
        BEH5p = MISSING_VALUE;
        MATS4s = MISSING_VALUE;
        BEL5e = MISSING_VALUE;
        MATS3s = MISSING_VALUE;
        X3v = MISSING_VALUE;
        Mw = MISSING_VALUE;
        CATS2D_3_DL = MISSING_VALUE;
        P_VSA_m_4 = MISSING_VALUE;
    }

    public void CalculateDescriptors(InsilicoMolecule Mol) {
        CalculateMw(Mol);
        CalculateMATS(Mol);
        CalculateALogP(Mol);

        DescriptorMLogP descriptorMLogP = new DescriptorMLogP();
        descriptorMLogP.CalculateAllDescriptors(Mol);

        MLogP = descriptorMLogP.MLogP;

        CalculateX3v(Mol);
        CalculateBEHBEL(Mol);
        CalculatePVSA(Mol);
        CalculateCats2D(Mol);
    }

    private void CalculatePVSA(InsilicoMolecule Mol) {

        P_VSA_m_4 = 0;
        short WEIGHT_M_IDX = 0;; short bin_size = 4;


        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            P_VSA_m_4 = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();


        // Calculate VSA
        double[] VSA = new double[nSK];

        for (int i=0; i<nSK; i++) {

            IAtom at = m.getAtom(i);

            // R - van der waals radius
            double vdwR = GetVdWRadius(m, at);
            if (vdwR == MISSING_VALUE) continue;

            double coef = 0;
            for (IAtom connAt : m.getConnectedAtomsList(at)) {

                double connAt_vdwR = GetVdWRadius(m, connAt);
                if (connAt_vdwR == MISSING_VALUE) continue;

                // refR - reference bond length
                double refR = this.GetRefBondLength(at, connAt);
                if (refR == MISSING_VALUE) continue;

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



        // Cycle for all found weighting schemes

        // Sets needed weights
        double[] w = Mass.getWeights(m);

        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE) continue;
            if ( bin_size == CalculateBin(w[i]))
                P_VSA_m_4 += VSA[i];
        }
    }

    private void CalculateX3v(InsilicoMolecule Mol) {

        int MaxPath = 3;
        X3v = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            X3v = -999;
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            X3v = -999;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);
        double[] curDescX = new double[MaxPath+1];
        double[] curDescXv = new double[MaxPath+1];
        double[] curDescXsol = new double[MaxPath+1];

        // checks for missing weights
        for (int qnumber : Qnumbers)
            if (qnumber == -999)
                return;
        for (double v : ValenceVD)
            if (v == -999)
                return;

        // clears VD matrix from linked F
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ((ConnAugMatrix[i][j]>0) && (ConnAugMatrix[j][j]==9))
                    VD[i]--;
            }


        for (int k=0; k<MaxPath; k++) {
            curDescX[k] = 0;
            curDescXv[k] = 0;
            curDescXsol[k] = 0;
        }

        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account

            // path 0
            curDescX[0] += Math.pow(VD[i], -0.5);
            curDescXv[0] += Math.pow(ValenceVD[i], -0.5);
            curDescXsol[0] += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);

            // path 1 - MaxPath
            for (int path=1; path<(MaxPath+1); path++) {

                if (curDescX[path] == -999) continue;

                IAtom at = m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (List<IAtom> curPath : CurPaths) {
//                    double prodX = 1;
                    double prodXv = 1;
//                    int prodQuantum = 1;
                    for (IAtom iAtom : curPath) {
                        int atIdx = m.indexOf(iAtom);
//                        if (ConnMatrix[atIdx][atIdx] == 9)
//                            continue; // F not taken into account
//                        prodX *= VD[atIdx];
                        prodXv *= ValenceVD[atIdx];
//                        prodQuantum *= Qnumbers[atIdx];
                    }
//                    curDescX[path] += Math.pow(prodX, -0.5);
                    curDescXv[path] += Math.pow(prodXv, -0.5);
//                    curDescXsol[path] += (1.00 / Math.pow(2.00, (double) (path + 1))) *
//                            ((double) prodQuantum) * Math.pow(prodX, -0.5);
                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath+1); i++) {
//            curDescX[i] /= 2;
            curDescXv[i] /= 2;
//            curDescXsol[i] /= 2;
        }

        X3v = curDescXv[3];

    }

    private void CalculateCats2D(InsilicoMolecule Mol) {

        CATS2D_3_DL = 0;
        int MAX_CATS_DISTANCE = 10;

        String TYPE_L = "L";
        String TYPE_D = "D";
        String[][] AtomCouples = {
                {TYPE_D, TYPE_L},
        };

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_3_DL = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DL = MISSING_VALUE;
            return;
        }
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DL = MISSING_VALUE;
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (String[] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[MAX_CATS_DISTANCE];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (isIn(CatsTypes[j], atomCouple[1])) {

                            if (TopoMat[i][j] < MAX_CATS_DISTANCE)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++) {
                if (i == 3)
                    CATS2D_3_DL = desc[i];
            }

        }


    }

    private void CalculateBEHBEL(InsilicoMolecule mol) {

        BEH5p = 0; BEL5e = 0;

        short WEIGHT_E_IDX = 2;
        short WEIGHT_P_IDX = 1;

        int MaxEig = 8;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH5p = MISSING_VALUE;
            BEL5e = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        BurdenMat = BurdenMatrix.getMatrix(m);

        int nSK = m.getAtomCount();


        // Sets needed weights
        double[] w = Polarizability.getWeights(m);

        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == MISSING_VALUE)
                MissingWeight = true;

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

        for (int i=0; i<MaxEig; i++) {
            double valH;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
            }

            if (i == 4)
                BEH5p = valH;

//            SetByName("BEH" + (i+1) + WEIGHT_SYMBOL[curWeight], valH);
//            SetByName("BEL" + (i+1) + WEIGHT_SYMBOL[curWeight], valL);
        }


        w = Electronegativity.getWeights(m);


        // If one or more weights are not available, sets all to missing value
        MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<MaxEig; i++) {
            double valL;
            if (i>(eigenvalues.length-1)) {
                valL = 0;
            } else {
                valL = eigenvalues[i];
            }
            if(i==4)
                BEL5e = valL;
        }
    }



    private void CalculateALogP(InsilicoMolecule Mol){
        ALogP = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ALogP = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);

        double[] Frags = acf.GetAllValues();

        for (double d : Frags)
            if (d == MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            ALogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }
    }

    private void CalculateMATS(InsilicoMolecule Mol) {

        MATS3s = 0;
        MATS4s = 0;

        int PARAMETER_LAG_03 = 3;
        int PARAMETER_LAG_04 = 4;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MATS3s = MISSING_VALUE;
            MATS4s = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS3s = MISSING_VALUE;
            MATS4s = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w;

        try {
            EState ES = new EState(Mol.GetStructure());
            w = ES.getIS();
        } catch (Exception e) {
            w = new double[nSK];
            for (int i=0; i<nSK; i++) w[i]= Descriptor.MISSING_VALUE;
        }

        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        // Calculates autocorrelations


        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {
            denom += Math.pow((w[i] - wA), 2);
            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_03) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }
        MATS3s = MoranAC;

        MoranAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {
            denom += Math.pow((w[i] - wA), 2);
            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_04) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }
        MATS4s = MoranAC;

    }

    private void CalculateMw(InsilicoMolecule Mol) {

        Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            Mw = -999;
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[] H = new int[nSK];

            for (int i=0; i<nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) {
                    log.warn(e.getMessage());
                }
            }

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");

            for (int i=0; i<nSK; i++) {
                if (wMass[i] == -999)
                    Mw = -999;
            }

            for (int i=0; i<nSK; i++) {
                if (Mw != -999) {
                    Mw += wMass[i];
                    if (H[i] > 0) {
                        Mw += HMass * H[i];
                    }
                }
            }

        } catch (Throwable e) {
            log.warn(e.getMessage());
            Mw = -999;
        }

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

            String TYPE_L = "L";
            String TYPE_D = "D";


            // Sets final types
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

    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }

    private int CalculateBin(double value) {

        // Ionization potential
        if (value < 1) return 1;
        if (value < 1.2) return 2;
        if (value < 1.6) return 3;
        if (value < 3) return 4;
        return 5;

    }



}
