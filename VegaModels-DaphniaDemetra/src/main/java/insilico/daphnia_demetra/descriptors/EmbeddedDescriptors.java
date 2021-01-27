package insilico.daphnia_demetra.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


import insilico.daphnia_demetra.descriptors.weights.BurdenMatrix;
import insilico.daphnia_demetra.descriptors.weights.Mass;
import insilico.daphnia_demetra.descriptors.weights.Polarizability;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    public double BEH1m;
    public double Eig1p;
    public double IC2;
    public double IDE;
    public double SRW5;
    public double T_F_CL;
    public double WA;
    public double Mw;

    public static final int MISSING_VALUE = -999;

//    Logger log = LoggerFactory.getLogger(EmbeddedDescriptors.class);

    private final static short WEIGHT_M_IDX = 0;
    private final static String[][] AtomCouples = {
            {"F", "Cl"}
    };


    public EmbeddedDescriptors() {
        BEH1m = -999;
        Eig1p = -999;
        IC2 = -999;
        IDE = -999;
        SRW5 = -999;
        T_F_CL = -999;
        WA = -999;
        Mw = -999;
    }

    public void CalculateAllEmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateSRW5(Mol);
        CalculateInformationContent(Mol);
        CalculateE1g1p(Mol);
        CalculateBEH1m(Mol);
        CalculateTFCL(Mol);
        CalculateWA(Mol);
        CalculateMw(Mol);
    }

    private void CalculateMw(InsilicoMolecule Mol){

        Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            WA = -999;
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
            WA = -999;
        }


    }

    private void CalculateWA(InsilicoMolecule Mol){

        WA = 0;

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            WA = -999;
            return;
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            WA = -999;
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            WeightsVertexDegree weightsVertexDegree = new WeightsVertexDegree();
            int[] VD = WeightsVertexDegree.getWeights(curMol, true);
//            int[] VD = VertexDegree.getWeights(curMol, true);

            // Wiener index
            double W = 0;
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    W += TopoMat[i][j];
                }

            WA = W/(nSK*(nSK-1));


        } catch (Throwable e) {
            log.warn(e.getMessage());
            WA = -999;
        }
    }

    private void CalculateTFCL(InsilicoMolecule mol) {

        T_F_CL = 0;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            T_F_CL = -999;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            T_F_CL = -999;
            return;
        }

        for (String[] atomCouple : AtomCouples) {
            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCouple[1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2)
                                T_F_CL += TopoMat[i][j];
                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (atomCouple[0].compareTo(atomCouple[1]) == 0) {
                T_F_CL /= 2;
            }
        }
    }

    private void CalculateSRW5(InsilicoMolecule mol) {

        int maxPath = 6;
        int SRWMaxPath =  maxPath - 1;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            SRW5 = -999;
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SRW5 = -999;
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int[] MolSRW = GetSRWToLag(SRWMaxPath, AdjMatDbl);
        SRW5 = MolSRW[SRWMaxPath];

    }


    private void CalculateInformationContent(InsilicoMolecule mol) {
        IC2 = 0;
        IDE = 0;
        int MaxPath = 2;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            IC2 = MISSING_VALUE;
            IDE = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] ConnMat;
        int[][] TopoDistMat;
        try {
            ConnMat = mol.GetMatrixConnectionAugmented();
            TopoDistMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            IC2 = MISSING_VALUE;
            IDE = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = WeightsVertexDegree.getWeights(m, false);
//        int[] VertexDeg = VertexDegree.getWeights(m, true);

        int[] TopoDistFreq = new int[nSK];  // frequencies of topological distances
        double TopoDistFreqSum = 0;
        for (int i=0; i<nSK; i++)
            TopoDistFreq[i] = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i+1; j<nSK; j++) {
                TopoDistFreq[TopoDistMat[i][j]]++;
                TopoDistFreqSum += TopoDistMat[i][j];
            }

        // IDE
        double denom = (double)nSK*(nSK-1)/2.00;
        for (int i=0; i<nSK; i++)
            if (TopoDistFreq[i]>0)
                IDE += ((double)TopoDistFreq[i]/denom) * (Math.log((double)TopoDistFreq[i]/denom));
        IDE = (-1.00 / Math.log(2)) * IDE;

        // Neighborhood indices

        double bic_denom = 0;
        for (IBond bnd : m.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);



            // Create belonging class for each atom(vertex)
        ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK);
        for (int i=0; i<nSK; i++) {
            IAtom atStart = m.getAtom(i);
            ArrayList<String> CurNeig = new ArrayList<>();
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if (TopoDistMat[i][j] == MaxPath) {
                    IAtom atEnd =  m.getAtom(j);
                    ShortestPaths sps = new ShortestPaths(m, atStart);
                    List<IAtom> sp = Arrays.asList(sps.atomsTo(atEnd));
                    // OLD DEPRECATED METHOD
//                        List<IAtom> sp = PathTools.getShortestPath(m, atStart, atEnd);

                    StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                    for (int k=0; k<(sp.size()-1); k++) {
                        int a = m.indexOf(sp.get(k));
                        int b = m.indexOf(sp.get(k + 1));
                        if (ConnMat[a][b] == 1)
                            bufPath.append("s");
                        if (ConnMat[a][b] == 2)
                            bufPath.append("d");
                        if (ConnMat[a][b] == 3)
                            bufPath.append("t");
                        if (ConnMat[a][b] == 1.5)
                            bufPath.append("a");
                        bufPath.append(sp.get(k + 1).getAtomicNumber());
                        bufPath.append("(").append(VertexDeg[m.indexOf(sp.get(k + 1))]).append(")");
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

        // Calculate IC

        for (int i=0; i<Gn.size(); i++)
            IC2 += ((double)Gn.get(i)/nSK) * (Math.log((double)Gn.get(i)/nSK));
        IC2 = (-1.00 / Math.log(2)) * IC2;

    }

    private void CalculateE1g1p(InsilicoMolecule mol) {
        Eig1p = 0;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            Eig1p = MISSING_VALUE;
            return;
        }

        double[][] ConnMatrix;
        try {
            ConnMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            Eig1p = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

//        WeightsPolarizability polarizability = new WeightsPolarizability();

        double[] w = Polarizability.getWeights(m);
//        double[] w = polarizability.getScaledWeights(m);
        double refW = Polarizability.GetPolarizability("C");
//        double refW = polarizability.getScaledWeight("C");

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

        Eig1p= eigenvalues[eigenvalues.length-1];

    }

    private void CalculateBEH1m(InsilicoMolecule mol) {
        BEH1m = 0;

        IAtomContainer m;

        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH1m = MISSING_VALUE;
            return;
        }

        double[][] BurdenMat;
        BurdenMat = BurdenMatrix.getMatrix(m);

        int nSK = m.getAtomCount();

//        ArrayList<Integer> weightList = BuildWeightListBurdenEigenvalue();

//        WeightsMass mass = new WeightsMass();
//        double[] w = mass.getScaledWeights(m);
        double[] w = Mass.getWeights(m);

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

        BEH1m = eigenvalues[eigenvalues.length-1];
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




}
