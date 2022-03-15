package insilico.aromatase_irfmn_tk.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.aromatase_irfmn_tk.descriptors.weights.ACF;
import insilico.aromatase_irfmn_tk.descriptors.weights.DescriptorMLogP;
import insilico.aromatase_irfmn_tk.descriptors.weights.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.aromatase_irfmn_tk.descriptors.weights.EStateCorrectForH;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    public double BEH2m;
    public double BEH2p;
    public double piPC08;
    public double piPC09;
    public double MLogP;
    public double H049;
    public double X5sol;
    public double nRNR2;
    public double ATSC4p;
    public double ATS5s;
    public double SsssN;
    public double SssNH;
    public double F02C_N;
    public double EEig3bo;

    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        MLogP = descriptorMLogP.MLogP;

        ACF acf  = new ACF();
        acf.Calculate(Mol);
        H049 = acf.GetByName("H-049").getValue();

        CalculateBEH(Mol);
        CalculateWAP(Mol);
        CalculateX5Sol(Mol);
        CalculatenRNR2(Mol);
        CalculateATS(Mol);
        CalculateEStates(Mol);
        CalculateF02C_N(Mol);
        CalculateEEig3bo(Mol);
    }

    private void CalculateEEig3bo(InsilicoMolecule Mol){
        EEig3bo = 0;
        int MinEig = 1;
        int MaxEig = 15;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            EEig3bo = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            EEig3bo = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            EEig3bo = MISSING_VALUE;
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
        for (int i=MinEig; i<=MaxEig; i++) {
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                if(i == 3)
                    EEig3bo = eigenvalues[idx];
        }

    }

    private void CalculateF02C_N(InsilicoMolecule Mol){

        F02C_N = 0;
        String[][] AtomCouples = {
                {"C", "N"}
        };

        short MAX_TOPO_DISTANCE = 2;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            F02C_N = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            F02C_N = MISSING_VALUE;
            return;
        }

        for (String[] atomCouple : AtomCouples) {

            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descF, 0);

            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCouple[1])) {

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descF[TopoMat[i][j] - 1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (atomCouple[0].compareTo(atomCouple[1]) == 0) {
                for (int i = 0; i < descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i = 0; i < descF.length; i++)
                if (i == 1)
                    F02C_N = descF[i];


        }
    }

    private void CalculateEStates(InsilicoMolecule Mol) {
        SsssN = 0;
        SssNH = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            SsssN = MISSING_VALUE;
            SssNH = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            log.warn(e.getMessage());
            SsssN = MISSING_VALUE;
            SssNH = MISSING_VALUE;
            return;
        }


        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("unable to get H count");
            }


            // Count bonds
            int nBnd=0, nSng = 0;
            for (IBond b : m.getConnectedBondsList(curAt)) {
                if (b.getFlag(CDKConstants.ISAROMATIC)) {
                    nBnd++;
                    continue;
                }
                if (b.getOrder() == IBond.Order.SINGLE) {
                    nSng++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.DOUBLE) {
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nBnd++;
                }
            }

            if (curAt.getSymbol().equalsIgnoreCase("N")) {


                if ((nBnd == 3) && (nSng == 3))
                    SsssN += es.getEState()[at];

                if ((nBnd == 2) && (nSng == 2) && (nH == 1))
                    SssNH += es.getEState()[at];
            }




        }


    }

    private void CalculateATS(InsilicoMolecule Mol){
        ATS5s = 0; ATSC4p = 0;
        short lag_p = 4;
        short lag_s = 5;

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            ATS5s = MISSING_VALUE; ATSC4p = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            ATS5s = MISSING_VALUE; ATSC4p = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w_s;
        double[] w_p;

        w_p = Polarizability.getWeights(m);

        try {
            EStateCorrectForH ES = new EStateCorrectForH(m);
            w_s = ES.getIS();

            // correction for compatibility with D7
            // H I-state is always 1
            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                    w_s[i] = 1;
            }
        } catch (Exception e) {
            w_s = new double[nSK];
            for (int i=0; i<nSK; i++)
                w_s[i]=Descriptor.MISSING_VALUE;
        }

        // Calculates weights averages
        double wA_P = 0;
        double wA_S = 0;
        for (int i=0; i<nSK; i++) {
            wA_P += w_p[i];
            wA_S += w_s[i];

        }
        wA_P = wA_P / ((double) nSK);
        wA_S = wA_S / ((double) nSK);

        // s
        double AC=0;
        double denom = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_s[i] - wA_S), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_s) {
                    AC += w_s[i] * w_s[j];
                }
        }


        // AC transformed in log form
        AC /= 2.0;
        ATS5s = Math.log(1 + AC);

        // p
        double ACS=0;
        denom = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_p[i] - wA_P), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_p) {
                    ACS += Math.abs((w_p[i]-wA_P) * (w_p[j]-wA_P));
                }
        }


        // AC transformed in log form
        ACS /= 2.0;
        ATSC4p = ACS;


    }

    private void CalculatenRNR2(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        nRNR2 = 0;
        DescriptorBlock descriptorBlock = new FunctionalGroups();
        descriptorBlock.Calculate(Mol);
        nRNR2 = descriptorBlock.GetByName("nRNR2").getValue();
    }

    private void CalculateBEH(InsilicoMolecule Mol){
        BEH2m = 0;
        BEH2p = 0;


        int MaxEig = 8;

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            log.warn(e.getMessage());
            BEH2m = MISSING_VALUE;
            BEH2p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat_m;
        double[][] BurdenMat_p;
        BurdenMat_m = BurdenMatrix.getMatrix(m);
        BurdenMat_p = BurdenMatrix.getMatrix(m);

        int nSK = m.getAtomCount();
        double[] w_p = Polarizability.getWeights(m);
        double[] w_m = Mass.getWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat_m[i][i] = w_m[i];
            BurdenMat_p[i][i] = w_p[i];
        }

        Matrix DataMatrix_m = new Matrix(BurdenMat_m);
        Matrix DataMatrix_p = new Matrix(BurdenMat_p);
        double[] eigenvalues_m;
        double[] eigenvalues_p;
        EigenvalueDecomposition ed_m = new EigenvalueDecomposition(DataMatrix_m);
        EigenvalueDecomposition ed_p = new EigenvalueDecomposition(DataMatrix_p);
        eigenvalues_m = ed_m.getRealEigenvalues();
        eigenvalues_p = ed_p.getRealEigenvalues();
        Arrays.sort(eigenvalues_m);
        Arrays.sort(eigenvalues_p);

        // p
        for (int i=0; i<MaxEig; i++) {
            if(i == 1) {
                double valH;
                if (i>(eigenvalues_p.length-1)) {
                    valH = 0;
                } else {
                    valH = eigenvalues_p[eigenvalues_p.length-1-i];
                }
                BEH2p = valH;
            }

        }

        // m
        for (int i=0; i<MaxEig; i++) {
            if(i == 1) {
                double valH;
                if (i>(eigenvalues_m.length-1)) {
                    valH = 0;
                } else {
                    valH = eigenvalues_m[eigenvalues_m.length-1-i];
                }
                BEH2m = valH;
            }

        }

    }


    private void CalculateX5Sol(InsilicoMolecule Mol) {
        X5sol = 0;
        short MaxPath = 5;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            X5sol = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            X5sol = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);

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

            curDescXsol[k] = 0;
        }

        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account

            // path 0

            curDescXsol[0] += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);

            // path 1 - MaxPath
            for (int path=1; path<(MaxPath+1); path++) {


                IAtom at = m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (List<IAtom> curPath : CurPaths) {
                    double prodX = 1;
                    double prodXv = 1;
                    int prodQuantum = 1;
                    for (IAtom iAtom : curPath) {
                        int atIdx = m.indexOf(iAtom);
                        prodX *= VD[atIdx];
                        prodXv *= ValenceVD[atIdx];
                        prodQuantum *= Qnumbers[atIdx];
                    }

                    curDescXsol[path] += (1.00 / Math.pow(2.00, (double) (path + 1))) *
                            ((double) prodQuantum) * Math.pow(prodX, -0.5);
                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath+1); i++) {
            curDescXsol[i] /= 2;
        }
        
        
        // Sets descriptors
        for (int i=0; i<(MaxPath+1); i++) {

            if(i == MaxPath)
                X5sol = curDescXsol[i];
        }
    }

    private void CalculateWAP(InsilicoMolecule Mol){
        piPC08 = 0;
        piPC09 = 0;
        int path8 = 8;
        int path9 = 9;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            piPC08 = MISSING_VALUE; piPC09 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            piPC08 = MISSING_VALUE; piPC09 = MISSING_VALUE;
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();

        //8
        double CurMPC=0;
        double[][] AtomPath = GetAtomsPaths(path8, m);

        for (int i=0; i<nSK; i++) {
            CurMPC += AtomPath[i][1];
        }

        CurMPC /= 2.0;
        CurMPC = Math.log(1+CurMPC);

        piPC08 = CurMPC;

        //9
        CurMPC=0;
        AtomPath = GetAtomsPaths(path9, m);

        for (int i=0; i<nSK; i++) {
            CurMPC += AtomPath[i][1];
        }

        CurMPC /= 2.0;
        CurMPC = Math.log(1+CurMPC);

        piPC09 = CurMPC;
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
