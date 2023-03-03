package insilico.carcinogenicity_sfi_classification.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.carcinogenicity_sfi_classification.descriptors.weights.*;
//import insilico.carcinogenicity_sfi_classification.descriptors.weights.*;
import insilico.core.descriptor.Descriptor;
//import insilico.core.descriptor.blocks.old.weight.Electronegativity;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;


public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    int MISSING_VALUE = -999;
    public double piPC10 = MISSING_VALUE;
    public double ATSC2p = MISSING_VALUE;
    public double EEig15bo = MISSING_VALUE;
    public double F1N_N= MISSING_VALUE;
    public double BEL3m = MISSING_VALUE;
    public double GATS1v = MISSING_VALUE;
    public double PW4 = MISSING_VALUE;
    public double MATS3s = MISSING_VALUE;
    public double BEH4e = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateWAP(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateBurdenEig(Mol);
        CalculateEAC(Mol);
        CalculateF1NN(Mol);
    }

    private void CalculateF1NN(InsilicoMolecule Mol){

        F1N_N = 0;
        String[][] AtomCouples = {
                {"N", "N"}
        };

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            F1N_N = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            F1N_N = MISSING_VALUE;
            return;
        }

        for (int d=0; d<AtomCouples.length; d++) {

            int[] descF = new int[10];
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouples[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouples[d][1])) {

                            if (TopoMat[i][j] <= 10) {
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouples[d][0].compareTo(AtomCouples[d][1]) == 0) {
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i=0; i<descF.length; i++)
                if(i == 0)
                    F1N_N = descF[i];
        }

    }

    private void CalculateEAC(InsilicoMolecule Mol){
        EEig15bo = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            EEig15bo = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            EEig15bo = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            EEig15bo = MISSING_VALUE;
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


        // EEig
        for (int i=1; i<=15; i++) {
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                if(i == 15)
                EEig15bo = eigenvalues[idx];
        }

    }

    private void CalculateBurdenEig(InsilicoMolecule Mol){
        BEH4e = 0; BEL3m = 0;
        short MaxEig = 8;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH4e = MISSING_VALUE; BEL3m = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        // Gets matrices
        double[][] BurdenMat_m;
        double[][] BurdenMat_e;
        try {
            BurdenMat_m = BurdenMatrix.getMatrix(m);
            BurdenMat_e = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            BEH4e = MISSING_VALUE; BEL3m = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w_e = Electronegativity.getWeights(m);
        double[] w_m = Mass.getWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat_m[i][i] = w_m[i];
            BurdenMat_e[i][i] = w_e[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix_m = new Matrix(BurdenMat_m);
        Matrix DataMatrix_e = new Matrix(BurdenMat_e);
        double[] eigenvalues_m;
        double[] eigenvalues_e;
        EigenvalueDecomposition ed_m = new EigenvalueDecomposition(DataMatrix_m);
        EigenvalueDecomposition ed_e = new EigenvalueDecomposition(DataMatrix_e);
        eigenvalues_m = ed_m.getRealEigenvalues();
        eigenvalues_e = ed_e.getRealEigenvalues();
        Arrays.sort(eigenvalues_m);
        Arrays.sort(eigenvalues_e);

        for (int i=0; i<MaxEig; i++) {
            double valL;
            if (i>(eigenvalues_m.length-1)) {
                valL = 0;
            } else {
                valL = eigenvalues_m[i];
            }
            if (i == 2)
                BEL3m = valL;

        }

        for (int i=0; i<MaxEig; i++) {
            double valH;
            if (i>(eigenvalues_e.length-1)) {
                valH = 0;
            } else {
                valH = eigenvalues_e[eigenvalues_e.length-1-i];
            }
            if (i == 3)
                BEH4e = valH;
        }

    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){
        ATSC2p = 0; GATS1v = 0; MATS3s = 0;
        double[] w_p, w_v, w_s;

        short lag_v = 1;
        short lag_p = 2;
        short lag_s = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ATSC2p = MISSING_VALUE; GATS1v = MISSING_VALUE; MATS3s = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            ATSC2p = MISSING_VALUE; GATS1v = MISSING_VALUE; MATS3s = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        w_p = Polarizability.getWeights(m);
        w_v = VanDerWaals.getWeights(m);;
        try {
            EState ES = new EState(Mol.GetStructure());
            w_s = ES.getIS();
        } catch (Exception e) {
            w_s = new double[nSK];
            for (int i=0; i<nSK; i++)
                w_s[i]=Descriptor.MISSING_VALUE;
        }

        double wA_p = 0;
        double wA_v = 0;
        double wA_s = 0;
        for (int i=0; i<nSK; i++) {
            wA_p += w_p[i];
            wA_v += w_v[i];
            wA_s += w_s[i];
        }

        wA_p = wA_p / ((double) nSK);
        wA_v = wA_v / ((double) nSK);
        wA_s = wA_s / ((double) nSK);

        double ACS=0;
        double denom = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_p[i] - wA_p), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_p) {
                    ACS += Math.abs((w_p[i]-wA_p) * (w_p[j]-wA_p));
                }
        }

        ACS /= 2.0;
        ATSC2p = ACS;

        double MoranAC=0;
        denom = 0;
        double delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_s[i] - wA_s), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_s) {
                    MoranAC += (w_s[i] - wA_s) * (w_s[j] - wA_s);
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

        double GearyAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_v[i] - wA_v), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_v) {
                    GearyAC += Math.pow((w_v[i] - w_v[j]), 2);
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

        GATS1v = GearyAC;

    }

    private void CalculateWAP(InsilicoMolecule Mol){
        PW4 = 0;
        piPC10 = 0;

        short pw_path = 4;
        short pi_path = 10;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            PW4 = MISSING_VALUE;
            piPC10 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            PW4 = MISSING_VALUE;
            piPC10 = MISSING_VALUE;
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();

        double CurPW=0;
        int[] AtomWalk = GetAtomsWalks(pw_path, AdjMatDbl);
        double[][] AtomPath = GetAtomsPaths(pw_path, m);

        for (int i=0; i<nSK; i++) {
            CurPW += (double)AtomPath[i][0] / (double)AtomWalk[i];
        }

        CurPW /= nSK;
        PW4 = CurPW;

        double CurMPC=0;
        AtomPath = GetAtomsPaths(pi_path, m);

        for (int i=0; i<nSK; i++) {
            CurMPC += AtomPath[i][1];
        }

        CurMPC /= 2.0;
        CurMPC = Math.log(1+CurMPC);
        piPC10 = CurMPC;

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
                for (List<IAtom> iAtoms : list) {
                    if (iAtoms.size() == (PathsOrder + 1)) {
                        PathNum++;

                        double curBO = 1;
                        for (int at_idx = 0; at_idx < iAtoms.size() - 1; at_idx++) {
                            IAtom at_1 = iAtoms.get(at_idx);
                            IAtom at_2 = iAtoms.get(at_idx + 1);
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
                        totBO += curBO;
                    }
                }
                paths[i][0] = PathNum;
                paths[i][1] = totBO;
            }
        }

        return paths;
    }
}
