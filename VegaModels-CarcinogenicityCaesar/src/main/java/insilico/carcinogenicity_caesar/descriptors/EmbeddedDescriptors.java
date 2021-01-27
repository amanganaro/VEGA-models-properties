package insilico.carcinogenicity_caesar.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.carcinogenicity_caesar.descriptors.weights.Polarizability;
import insilico.carcinogenicity_caesar.descriptors.weights.VertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;

import java.util.Arrays;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    public double PW5;
    public double DDr6;
    public double MATS2p;
    public double EEig10ed;
    public double ESpm11ed;
    public double ESpm9dm;
    public double GGI2;
    public double JGI6;

    public EmbeddedDescriptors(){

        DDr6 = MISSING_VALUE;
        MATS2p = MISSING_VALUE;
        EEig10ed = MISSING_VALUE;
        ESpm11ed = MISSING_VALUE;
        ESpm9dm = MISSING_VALUE;
        GGI2 = MISSING_VALUE;
        JGI6 = MISSING_VALUE;
        PW5 = MISSING_VALUE;
    }

    public void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateDDr6(Mol);
        CalculatePW5(Mol);
        CalculateEdgeAdjacency(Mol);
        CalculateTopologicalCharge(Mol);
        CalculateMATS2p(Mol);
    }

    private void CalculateTopologicalCharge(InsilicoMolecule mol) {
        GGI2 = 0;
        JGI6 = 0;
        int MaxPath = 8;

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            GGI2 = MISSING_VALUE;
            JGI6 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMat;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            GGI2 = MISSING_VALUE;
            JGI6 = MISSING_VALUE;
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[][] AdjMat = AdjacencyMatrix.getMatrix(curMol);
            int[] VertexDegrees = VertexDegree.getWeights(curMol, false);

            Matrix mAdj = new Matrix(AdjMat.length, AdjMat[0].length);
            for (int i=0; i<AdjMat.length; i++)
                for (int j=0; j<AdjMat[0].length; j++)
                    mAdj.set(i, j, (double)AdjMat[i][j]);

            Matrix mRecSqrDist = new Matrix(TopoMat.length, TopoMat[0].length);
            for (int i=0; i<TopoMat.length; i++)
                for (int j=0; j<TopoMat[0].length; j++) {
                    double CurVal = 0;
                    if (TopoMat[i][j] != 0)
                        CurVal = 1 / Math.pow((double)TopoMat[i][j], 2);
                    mRecSqrDist.set(i, j, CurVal);
                }

            Matrix multMat = mAdj.times(mRecSqrDist);

            double[][] CTMatrix = new double[nSK][nSK];
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++) {
                    if (i == j)
                        CTMatrix[i][j] = VertexDegrees[i];
                    else {
                        CTMatrix[i][j] = multMat.get(i, j) - multMat.get(j, i);
                    }
                }

            int[] PathCount = new int[MaxPath];
            double[] GGIval = new double[MaxPath];
            for (int i=0; i<MaxPath; i++) {
                PathCount[i] = 0;
                GGIval[i] = 0;
            }

            for (int i=0; i<nSK; i++)
                for (int j=i; j<nSK; j++) {

                    int CurPath = TopoMat[i][j];
                    if  ((CurPath>0) && (CurPath<=MaxPath)) {
                        PathCount[CurPath-1]++;
                        GGIval[CurPath-1] += Math.abs(CTMatrix[i][j]);
                    }
                }

            // Sets descriptors
//            double[] GGI = new double[MaxPath + 1];
//            double[] JGI = new double[MaxPath + 1];

            if (GGIval[1] > 0)
                GGI2 = GGIval[1];
            else GGI2 = 0;
            if(PathCount[5] != 0)
                JGI6 = GGIval[5] / (double) PathCount[5];
            else JGI6 = 0;

        } catch (Throwable e) {
            GGI2 = MISSING_VALUE;
            JGI6 = MISSING_VALUE;
        }
    }

    private void CalculateEdgeAdjacency(InsilicoMolecule mol) {
        // EEig + 10 + ed
        // EEig + 11 + ed
        // ESpm + 9 + dm
        EEig10ed = 0; ESpm11ed = 0; ESpm9dm = 0;

        int MinEig = 1;
        int MaxEig = 15;
        int MinSM = 1;
        int MaxSM = 15;

        short WEIGHT_X_IDX = 0;
        short WEIGHT_D_IDX = 1;

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            EEig10ed = MISSING_VALUE; ESpm11ed = MISSING_VALUE; ESpm9dm = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            EEig10ed = MISSING_VALUE; ESpm11ed = MISSING_VALUE; ESpm9dm = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            EEig10ed = MISSING_VALUE; ESpm11ed = MISSING_VALUE; ESpm9dm = MISSING_VALUE;
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeDegreeMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                EdgeDegreeMat[i][j] = EdgeAdjMat[i][j][1];

        DataMatrix = new Matrix(EdgeDegreeMat);

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        double[] EEig = new double[MaxEig + 1];

        // EEig
        for (int i=MinEig; i<=MaxEig; i++) {
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                EEig[i] = eigenvalues[idx];
            else
                EEig[i] = 0;
        }

        EEig10ed = EEig[10];

//        ESpm11ed
        double[] Espm = new double[MaxSM + 1];

        for (int i=MinSM; i<=MaxSM; i++) {
            double curSM = 0;
            for (int k=(eigenvalues.length-1); k>=0; k--) {
                curSM += Math.pow(eigenvalues[k], (i));
            }
            curSM = Math.log(1 + curSM);

            Espm[i] = curSM;
        }

        ESpm11ed = Espm[11];

        DataMatrix = null;
        double[][] EdgeDipoleMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                EdgeDipoleMat[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<m.getBondCount(); i++) {
            IAtom a =  m.getBond(i).getAtom(0);
            IAtom b =  m.getBond(i).getAtom(1);

            double CurVal = GetDipoleMoment(m, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(m, b, a);
            EdgeDipoleMat[i][i] = CurVal;
        }

        DataMatrix = new Matrix(EdgeDipoleMat);

        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        Espm = new double[MaxSM+1];

        // Spectral moment
        for (int i=MinSM; i<=MaxSM; i++) {
            double curSM = 0;
            for (int k=(eigenvalues.length-1); k>=0; k--) {
                curSM += Math.pow(eigenvalues[k], (i));
            }
            curSM = Math.log(1 + curSM);

            Espm[i] = curSM;
        }

        ESpm9dm = Espm[9];
    }

    private void CalculateMATS2p(InsilicoMolecule Mol){
        MATS2p = 0;
        int PARAMETER_LAG_02 = 2;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MATS2p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS2p = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = Polarizability.getWeights(m);

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
                if (TopoMatrix[i][j] == PARAMETER_LAG_02) {
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
        MATS2p = MoranAC;
    }

    private void CalculateDDr6(InsilicoMolecule Mol) {

        int MinRingSize = 3;
        int MaxRingSize = 12;
        DDr6 = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            DDr6 = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        double[][] DDrMatrix = null;
        try {
            DDrMatrix = Mol.GetMatrixDistanceDetour();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            DDr6 = MISSING_VALUE;
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
            DDr6 = DDr[6];

        } catch (Exception e) {
            DDr6 = MISSING_VALUE;
        }
    }

    private void CalculatePW5(InsilicoMolecule Mol) {

        PW5 = 0;

        int PATH_PW = 5;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            PW5 = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            PW5 = MISSING_VALUE;
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

            PW5 = CurPW / nSK;
        }
    }




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


    private double GetDipoleMoment(IAtomContainer CurMol, IAtom at1, IAtom at2) {

        String a = at1.getSymbol();
        String b = at2.getSymbol();

        // C - something
        if (a.equalsIgnoreCase("C")) {

            // C-F
            if (b.equalsIgnoreCase("F")) {
                return 1.51;
            }

            // C-Cl , C(Cl)-Cl , C(Cl)(Cl)-Cl
            if (b.equalsIgnoreCase("Cl")) {
                int nCl=0;
                for (IAtom at : CurMol.getConnectedAtomsList(at1)) {
                    if (at.getSymbol().equalsIgnoreCase("Cl"))
                        nCl++;
                }
                if (nCl==1)
                    return 1.56;
                if (nCl==2)
                    return 1.20;
                if (nCl==3)
                    return 0.83;
            }

            // C-Br
            if (b.equalsIgnoreCase("Br")) {
                return 1.48;
            }

            // C-I
            if (b.equalsIgnoreCase("I")) {
                return 1.29;
            }

            // C-N , C=N , C#N
            if (b.equalsIgnoreCase("N")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.4;
                if (ord == IBond.Order.DOUBLE)
                    return 0.9;
                if (ord == IBond.Order.TRIPLE)
                    return 3.6;
            }

            // C-O , C=O
            if (b.equalsIgnoreCase("O")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.86;
                if (ord == IBond.Order.DOUBLE)
                    return 2.4;
            }

            // C-S , C=S
            if (b.equalsIgnoreCase("S")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 2.95;
                if (ord == IBond.Order.DOUBLE)
                    return 2.8;
            }

        }


        // N-O , N-[O-] , N=O
        if ((a.equalsIgnoreCase("N")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nH=0;
            try { nH = at2.getImplicitHydrogenCount(); } catch (Exception e) {}
            int nConn = CurMol.getConnectedBondsCount(at2) + nH;
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 3.2;
            if ((ord == IBond.Order.SINGLE) && (nConn==1))
                return 2.0;
//                return 0.3;
            if ((ord == IBond.Order.DOUBLE) && (nConn==1))
                return 2.0;
        }


        // S-[O-]
        if ((a.equalsIgnoreCase("S")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nConn = CurMol.getConnectedBondsCount(at2);
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 2.9;
        }


        // C(*)(*)-C(*)(*)(*) , C(*)(*)-C , CC(*)(*)(*)
        if ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("C"))) {

            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;

            int nH1=0, nH2=0;
            try {
                nH1 = at1.getImplicitHydrogenCount();
            } catch (Exception E) {}
            try {
                nH2 = at2.getImplicitHydrogenCount();
            } catch (Exception E) {}

            int nConn1 = CurMol.getConnectedBondsCount(at1) + nH1;
            int nConn2 = CurMol.getConnectedBondsCount(at2) + nH2;

            if ((nConn1==3) && (nConn2==4))
                return 0.68;
            if ((nConn1==3) && (nConn2==2))
                return 1.15;
            if ((nConn1==2) && (nConn2==4))
                return 1.48;
        }

        return 0;
    }



}
