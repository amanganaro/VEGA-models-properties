package insilico.watersolubility.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAMeylanLogPAdditionalFragments;
import insilico.core.alerts.builders.SAMeylanLogPCorrectionFragments;
import insilico.core.alerts.builders.SAMeylanLogPFragments;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
//import insilico.core.wa.VertexDegree;
import insilico.watersolubility.descriptors.weights.*;
import insilico.watersolubility.descriptors.weights.VertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.watersolubility.descriptors.weights.*;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;

import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.*;

import java.util.*;

@Slf4j
public class EmbeddedDescriptors {

    public static final int MISSING_VALUE = -999;

    public double ALogP;
    public double piPC8;
    public double SEigm;
    public double MATS1v;
    public double PW6;
    public double P_VSA_I_4;
    public double LogP;
    public double nR10;
    public double CATS2D_4_DL;
    public double BEH2p;
    public double BIC3;
    public double GATS2m;
    public double H050;
    public double DDr3;
    public double CATS2D_7_DL;

    private double[][] ConnAugMatrix;

    public EmbeddedDescriptors() {
        ALogP = MISSING_VALUE;
        piPC8 = MISSING_VALUE;
        PW6 = MISSING_VALUE;
        SEigm = MISSING_VALUE;
        MATS1v = MISSING_VALUE;
        GATS2m = MISSING_VALUE;
        P_VSA_I_4 = MISSING_VALUE;
        LogP = MISSING_VALUE;
        nR10 = MISSING_VALUE;
        CATS2D_4_DL = MISSING_VALUE;
        CATS2D_7_DL = MISSING_VALUE;
        BEH2p = MISSING_VALUE;
        BIC3 = MISSING_VALUE;
        H050 = MISSING_VALUE;
        DDr3 = MISSING_VALUE;
    }


    public void CalculateAllEmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateALogP(Mol);
        CalculateWalkAndPath(Mol);
        CalculateSEigm(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateP_VSA_1_4(Mol);
        CalculateMeylanLogP(Mol);
        CalculatenR10(Mol);
        CalculateCATS2D(Mol);
        CalculateBEH2p(Mol);
        CalculateBIC3(Mol);
        CalculateH050(Mol);
        CalculateDDr3(Mol);
    }

    private void CalculateDDr3(InsilicoMolecule Mol) {

        int MinRingSize = 3;
        int MaxRingSize = 12;
        DDr3 = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            DDr3 = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        double[][] DDrMatrix = null;
        try {
            DDrMatrix = Mol.GetMatrixDistanceDetour();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            DDr3 = MISSING_VALUE;
            return;
        }

        double[] DDr = new double[MaxRingSize+1];

        try {
            RingSet MolRings = Mol.GetAllRings();

            Arrays.fill(DDr, 0);

            for (int i=0; i< MolRings.getAtomContainerCount(); i++) {
                IRing r = (IRing) MolRings.getAtomContainer(i);
                int rSize = r.getAtomCount();

                if ((rSize >= MinRingSize) && (rSize <= MaxRingSize)) {
                    for (IAtom at : r.atoms()) {
                        int atNum = m.indexOf(at);
                        double rowSum = 0;
                        for (int k=0; k<DDrMatrix[atNum].length; k++)
                            rowSum += DDrMatrix[atNum][k];
                        DDr[rSize]  += rowSum;
                    }
                }

            }
            DDr3 = DDr[3];

        } catch (Exception e) {
            DDr3 = MISSING_VALUE;
        }
    }

    private void CalculateH050(InsilicoMolecule Mol) {
        H050 = 0;
        try {
            ACF acf = new ACF();
            acf.Calculate(Mol);
            H050 = acf.GetByName("H-050").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
            H050 = MISSING_VALUE;
        }


    }

    private void CalculateBIC3(InsilicoMolecule Mol) {

        BIC3 = 0;
        int MaxLag = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BIC3 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] ConnMat;
        int[][] TopoDistMat;
        try {
            ConnMat = Mol.GetMatrixConnectionAugmented();
            TopoDistMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            BIC3 = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);

        // Neighborhood indices

        double bic_denom = 0;
        for (IBond bnd : m.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag=1; CurLag<=MaxLag; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK);
            for (int i=0; i<nSK; i++) {
                IAtom atStart = m.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (TopoDistMat[i][j] == CurLag) {
                        IAtom atEnd =  m.getAtom(j);
//                        ShortestPaths sps = new ShortestPaths(m, atStart);
//                        List<IAtom> sp = Arrays.asList(sps.atomsTo(atEnd));
                        // OLD DEPRECATED METHOD
                        List<IAtom> sp = getShortestPath(m, atStart, atEnd);

                        String bufPath = "" + sp.get(0).getAtomicNumber();
                        for (int k=0; k<(sp.size()-1); k++) {
                            int a = m.indexOf(sp.get(k));
                            int b = m.indexOf(sp.get(k + 1));
                            if (ConnMat[a][b] == 1)
                                bufPath += "s";
                            if (ConnMat[a][b] == 2)
                                bufPath += "d";
                            if (ConnMat[a][b] == 3)
                                bufPath += "t";
                            if (ConnMat[a][b] == 1.5)
                                bufPath += "a";
                            bufPath += sp.get(k+1).getAtomicNumber();
                            bufPath += "(" + VertexDeg[m.indexOf(sp.get(k + 1))] + ")";
                        }
                        CurNeig.add(bufPath);
                    }
                }
                Collections.sort(CurNeig);
                NeigList.add(CurNeig);
            }

            // Calculates equivalence classes
            ArrayList<ArrayList<String>> G = new ArrayList<>();
            ArrayList<Integer> Gn = new ArrayList<>();
            for (int i=0; i<nSK; i++) {
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
                ic += ((double)Gn.get(i)/nSK) * (Math.log((double)Gn.get(i)/nSK));
            ic = (-1.00 / Math.log(2)) * ic;


            BIC3 = bic_denom==0 ? 0 : ic / bic_denom;


        }
    }

    private void CalculateBEH2p(InsilicoMolecule Mol) {

        BEH2p = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH2p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int nSK = m.getAtomCount();
        double[][] BurdenMat = new double[nSK][nSK];
        double[][] ConnMat = ConnectionAugMatrix.getMatrix(m);

        for(int i = 0; i < nSK; ++i) {
            for(int j = 0; j < nSK; ++j) {
                if (i == j) {
                    BurdenMat[i][j] = ConnMat[i][j];
                } else if (ConnMat[i][j] > 0.0D) {
                    BurdenMat[i][j] = Math.sqrt(ConnMat[i][j]);
//                    }
                } else {
                    BurdenMat[i][j] = 0.001D;
                }
            }
        }

        // Sets needed weights
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

        double valH, valL;
        if (1>(eigenvalues.length-1)) {
            valH = 0;
        } else {
            valH = eigenvalues[eigenvalues.length-2];
        }

        BEH2p = valH;

    }

    private void CalculateCATS2D(InsilicoMolecule Mol) {


        CATS2D_4_DL = 0;
        CATS2D_7_DL = 0;
        int MAX_CATS_DISTANCE = 10;

        String TYPE_A = "A";
        String TYPE_P = "P";
        String TYPE_N = "N";
        String TYPE_L = "L";
        String TYPE_D = "D";
        String[][] AtomCouples = {
                {TYPE_D, TYPE_L},
        };

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_4_DL = MISSING_VALUE;
            CATS2D_7_DL = MISSING_VALUE;
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
            CATS2D_4_DL = MISSING_VALUE;
            CATS2D_7_DL = MISSING_VALUE;
            return;
        }
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_4_DL = MISSING_VALUE;
            CATS2D_7_DL = MISSING_VALUE;
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (int d=0; d < AtomCouples.length; d++) {

            int descT = 0;
            int[] desc = new int[MAX_CATS_DISTANCE];
            Arrays.fill(desc, 0);

            for (int i=0; i<nSK; i++) {
                if ( isIn(CatsTypes[i], AtomCouples[d][0]) ) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if ( isIn(CatsTypes[j], AtomCouples[d][1]) ) {

                            if (TopoMat[i][j] < MAX_CATS_DISTANCE)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i=0; i<desc.length; i++) {
                if (i == 4)
                    CATS2D_4_DL = desc[i];
                if (i == 7)
                    CATS2D_7_DL = desc[i];
            }

        }



    }


    public void CalculatenR10(InsilicoMolecule Mol) {
        nR10 = 0;
        int RING_SIZE = 10;
        int RING_SIZE_MIN = 3;
        int RING_SIZE_MAX = 12;

        try {

            int nSizes = RING_SIZE_MAX - RING_SIZE_MIN + 1;
            int[] RingCount = new int[nSizes];
            int[] RingSize = new int[nSizes];
            for (int i=0; i<nSizes; i++) {
                RingSize[i] = RING_SIZE_MIN + i;
                RingCount[i] = 0;
            }

            IRingSet allRings = Mol.GetAllRings();
            Iterator<IAtomContainer> RingsIterator = allRings.atomContainers().iterator();
            while (RingsIterator.hasNext()) {
                IRing ring = (IRing)RingsIterator.next();
                for (int i=0; i<nSizes; i++) {
                    if (ring.getAtomCount() == RingSize[i])
                        RingCount[i]++;
                }
            }

            nR10 = RingCount[7];

//            for (int i=0; i<nSizes; i++)
//                SetByName("nR" + RingSize[i], RingCount[i]);

        } catch (Throwable e) {
            nR10 = MISSING_VALUE;
        }

    }

    private void CalculateMeylanLogP(InsilicoMolecule mol) {
        LogP = 0;

        try {
            IAtomContainer m = mol.GetStructure();

            SAMeylanLogPFragments Frags = new SAMeylanLogPFragments();
            AlertList f_base = Frags.Calculate(mol);
            SAMeylanLogPAdditionalFragments FragsAdd = new SAMeylanLogPAdditionalFragments();
            AlertList f_add = FragsAdd.Calculate(mol);
            SAMeylanLogPCorrectionFragments FragsCorr = new SAMeylanLogPCorrectionFragments();
            AlertList f_corr = FragsCorr.Calculate(mol);

            double Coefficient = Frags.GetCoefficient() + FragsAdd.GetCoefficient();
            double Correction = FragsCorr.GetCoefficient();

            // main eq
            LogP = 0.2290 + Coefficient + Correction;

            // lower bound
            if (LogP < -5.0)
                LogP = -5.0;

        } catch (InitFailureException | InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
        }
    }

    private void CalculateP_VSA_1_4(InsilicoMolecule Mol) {
        P_VSA_I_4 = 0;
        short WEIGHT_I_IDX = 2; short bin_size = 4;




        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            P_VSA_I_4 = MISSING_VALUE;
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
        double[] w = IonizationPotential.getWeightsNormalized(m);

//        int bins = BuildBinSize(curWeight);


        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE) continue;
            if ( bin_size == CalculateBin(w[i]))
                P_VSA_I_4 += VSA[i];
        }

    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){
        MATS1v = 0;
        GATS2m = 0;
        String PARAMETER_WEIGHT_M = "weightm";
        String PARAMETER_WEIGHT_P = "weightp";
        int PARAMETER_LAG_01 = 1;
        int PARAMETER_LAG_02 = 2;
        short WEIGHT_M_IDX = 0;
        short WEIGHT_V_IDX = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MATS1v = MISSING_VALUE;
            GATS2m = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS1v = MISSING_VALUE;
            GATS2m = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = VanDerWaals.getWeights(m);

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
                if (TopoMatrix[i][j] == PARAMETER_LAG_01) {
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
        MATS1v = MoranAC;

        w = Mass.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        // Calculates autocorrelations

            double GearyAC=0;
            denom = 0;
            delta = 0;

            for (int i=0; i<nSK; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == PARAMETER_LAG_02) {
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

            // Sets descriptor
           GATS2m = GearyAC;
    }

    private void CalculateSEigm(InsilicoMolecule Mol){
        SEigm = 0;
        short WEIGHT_M_IDX = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SEigm = MISSING_VALUE;
            return;
        }

        double[][] ConnMatrix;
        try {
            ConnMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SEigm = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = Mass.getWeights(m);
        double refW = Mass.GetMass("C");

        double[][] EigMat = new double[nSK][nSK];
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {

                if (i==j) {

                    EigMat[i][j] = 1 - (refW / w[i]);

                } else {

                    // builds shortest path between i and j
                    IAtom at1 = m.getAtom(i);
                    IAtom at2 =  m.getAtom(j);
                    ShortestPaths sp = new ShortestPaths(m, at1);
                    List<IAtom> Path = Arrays.asList(sp.atomsTo(at2));

                    // OLD DEPRECATED METHOD
//                        List<IAtom> Path = PathTools.getShortestPath(m, at1, at2);

                    double val = 0;
                    for (int k=0; k<(Path.size()-1); k++) {
                        int a1 = m.indexOf(Path.get(k));
                        int a2 = m.indexOf(Path.get(k + 1));
                        double bond = ConnMatrix[a1][a2];
                        val += (1 / bond) * (Math.pow(refW, 2) / (w[a1] * w[a2]) );
                    }

                    EigMat[i][j] = val;

                }
            }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(EigMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);


        for (double eigenvalue : eigenvalues) {
            SEigm += eigenvalue;
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

    private void CalculateWalkAndPath(InsilicoMolecule Mol) {

        piPC8 = 0;
        PW6 = 0;

        int PATH_PW = 6;
        int PATH_PI = 8;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            piPC8 = MISSING_VALUE; PW6 = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            piPC8 = MISSING_VALUE; PW6 = MISSING_VALUE;
            return;
        }

        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();

        // PW6
        for (int curPath = 1; curPath<= PATH_PW; curPath++) {

            double CurPW = 0;
            int[] AtomWalk = GetAtomsWalks(curPath, AdjMatDbl);
            double[][] AtomPath = GetAtomsPaths(curPath, m);

            for (int i=0; i<nSK; i++) {
                CurPW += AtomPath[i][0] / (double)AtomWalk[i];
            }

            PW6 = CurPW / nSK;
        }

        // piPC8
        for (int curPath = 1; curPath<= PATH_PI; curPath++) {

            double CurMPC=0;
            double[][] AtomPath = GetAtomsPaths(curPath, m);

            for (int i=0; i<nSK; i++) {
                CurMPC += AtomPath[i][1];
            }

            CurMPC /= 2.0;
            piPC8 = Math.log(1+CurMPC);
        }
    }


    private int[] GetAtomsWalks(int WalksOrder, double[][] AdjMatrix) {

        int[] walks = new int[AdjMatrix.length];

        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);

        for (int k=1; k<WalksOrder; k++) {
            mWalks = mWalks.times(mAdj);
        }

        for (int i=0; i<AdjMatrix.length; i++) {
            int CurSum = 0;
            for (int j=0; j<AdjMatrix[0].length; j++)
                CurSum += mWalks.get(i, j);
            walks[i] = CurSum;
        }

        return walks;
    }

    // index 0 = count of atom paths
    // index 1 = sum of paths weighted by bond order
    private double[][] GetAtomsPaths(int PathsOrder, IAtomContainer Mol) {

        int nSK = Mol.getAtomCount();
        double[][] paths = new double[nSK][2];

        for (int i=0; i<nSK; i++) {
            IAtom a = Mol.getAtom(i);
            int PathNum = 0;
            double totBO = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j)
                    continue;
                IAtom b =  Mol.getAtom(j);
                List<List<IAtom>> list= PathTools.getAllPaths(Mol, a, b);
                for (int k=0; k<list.size(); k++) {
                    if (list.get(k).size() == (PathsOrder+1)) {
                        PathNum++;

                        double curBO=1;
                        for (int at_idx=0; at_idx<list.get(k).size()-1; at_idx++) {
                            IAtom at_1 = list.get(k).get(at_idx);
                            IAtom at_2 = list.get(k).get(at_idx+1);
                            IBond curBond = Mol.getBond(at_2, at_1);
                            if (curBond.getFlag(CDKConstants.ISAROMATIC))
                                curBO *= 1.5;
                            else {
                                if (curBond.getOrder() == IBond.Order.SINGLE) curBO *= 1;
                                if (curBond.getOrder() == IBond.Order.DOUBLE) curBO *= 2;
                                if (curBond.getOrder() == IBond.Order.TRIPLE) curBO *= 3;
                                if (curBond.getOrder() == IBond.Order.QUADRUPLE) curBO *= 4;
                            }
                        }
                        totBO+=curBO;
                    }
                }
                paths[i][0] = PathNum;
                paths[i][1] = totBO;
            }
        }

        return paths;
    }


    private int[] GetSRWToLag(int Lag, double[][] AdjMatrix) {

        int nSK = AdjMatrix.length;
        int[] MolSRW = new int[Lag+1];

        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);

        for (int k=1; k<(Lag+1); k++) {

            MolSRW[k] = 0;
            for (int i=0; i<nSK; i++) {
                MolSRW[k] += mWalks.get(i, i);
            }

            mWalks = mWalks.times(mAdj);
        }

        return MolSRW;
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

    private int CalculateBin(double value) {

        // Ionization potential
        if (value < 1) return 1;
        if (value < 1.15) return 2;
        if (value < 1.25) return 3;
        return 4;


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

    private static List<IAtom> getShortestPath(IAtomContainer atomContainer, IAtom start, IAtom end) {
        int natom = atomContainer.getAtomCount();
        int endNumber = atomContainer.indexOf(end);
//        int endNumber = atomContainer.getAtomNumber(end);
        int startNumber = atomContainer.indexOf(start);
//        int startNumber = atomContainer.getAtomNumber(start);
        int[] dist = new int[natom];
        int[] previous = new int[natom];
        for (int i = 0; i < natom; i++) {
            dist[i] = 99999999;
            previous[i] = -1;
        }
        dist[atomContainer.indexOf(start)] = 0;
//        dist[atomContainer.getAtomNumber(start)] = 0;

        List<IAtom> Slist = new ArrayList<IAtom>();
        List<Integer> Qlist = new ArrayList<Integer>();
        for (int i = 0; i < natom; i++) Qlist.add(i);

        while (true) {
            if (Qlist.size() == 0) break;

            // extract min
            int u = 999999;
            int index = 0;
            for (Integer tmp : Qlist) {
                if (dist[tmp] < u) {
                    u = dist[tmp];
                    index = tmp;
                }
            }
            Qlist.remove(Qlist.indexOf(index));
            Slist.add(atomContainer.getAtom(index));
            if (index == endNumber) break;

            // relaxation
            List<IAtom> connected = atomContainer.getConnectedAtomsList(atomContainer.getAtom(index));
            for (IAtom aConnected : connected) {
                int anum = atomContainer.indexOf(aConnected);
//                int anum = atomContainer.getAtomNumber(aConnected);
                if (dist[anum] > dist[index] + 1) { // all edges have equals weights
                    dist[anum] = dist[index] + 1;
                    previous[anum] = index;
                }
            }
        }

        ArrayList<IAtom> tmp = new ArrayList<IAtom>();
        int tmpSerial = endNumber;
        while (true) {
            tmp.add(0, atomContainer.getAtom(tmpSerial));
            tmpSerial = previous[tmpSerial];
            if (tmpSerial == startNumber) {
                tmp.add(0, atomContainer.getAtom(tmpSerial));
                break;
            }
        }
        return tmp;
    }


}
