package insilico.daphnia_noec.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.weight.VertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.daphnia_noec.descriptors.weights.*;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;

@Slf4j
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    public double Mw = MISSING_VALUE;
    public double OPerc = MISSING_VALUE;
    public double MLogP = MISSING_VALUE;
    public double nArNH2 = MISSING_VALUE;
    public double nS = MISSING_VALUE;
    public double MATS5s = MISSING_VALUE;
    public double MATS6m = MISSING_VALUE;
    public double EEig7dm = MISSING_VALUE;
    public double CATS2D_7_DL = MISSING_VALUE;
    public double C_026 = MISSING_VALUE;
    public double JGI3 = MISSING_VALUE;

    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";

    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };
    private double[][] ConnAugMatrixCATS;


    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    public void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateConstitutional(Mol);
        CalculateMLogP(Mol);
        CalculateFG(Mol);
        CalculateACF(Mol);
        CalculateTopologicalCharge(Mol);
        CalculateCATS2(Mol);
        CalculateEAC(Mol);
        CalculateAutocorrelation(Mol);
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol){
        MATS5s = 0; MATS6m = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            MATS5s = MISSING_VALUE; MATS6m = MISSING_VALUE;
            return;
        }

        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS5s = MISSING_VALUE; MATS6m = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        double[] w = Mass.getWeights(m);
        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == MISSING_VALUE)
                MissingWeight = true;

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 6) {
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

        MATS6m = MoranAC;

        try {
            EState ES = new EState(Mol.GetStructure());
            w = ES.getIS();
        } catch (Exception e) {
            w = new double[nSK];
            for (int i=0; i<nSK; i++) w[i]=MISSING_VALUE;
        }
        for (int i=0; i<nSK; i++)
            if (w[i] == MISSING_VALUE)
                MissingWeight = true;

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 5) {
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

        MATS5s = MoranAC;

    }


    private void CalculateEAC(InsilicoMolecule Mol){
        EEig7dm = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            EEig7dm = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            EEig7dm = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            EEig7dm = MISSING_VALUE;
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeDipoleMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                EdgeDipoleMat[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<m.getBondCount(); i++) {
            for (int j=0; j<m.getBondCount(); j++) {
                if (EdgeDipoleMat[i][j] != 0) {
                    IAtom a =  m.getBond(i).getAtom(0);
                    IAtom b =  m.getBond(i).getAtom(1);
                    double CurVal = GetDipoleMoment(m, a, b);
                    if (CurVal == 0)
                        CurVal = GetDipoleMoment(m, b, a);
                    EdgeDipoleMat[i][j] = CurVal;
                }
            }
        }

        DataMatrix = new Matrix(EdgeDipoleMat);

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);


        // EEig
        for (int i=1; i<=15; i++) {
            if(i == 7) {
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    EEig7dm = eigenvalues[idx];
                else
                    EEig7dm = 0;
            }
        }
    }

    private void CalculateCATS2(InsilicoMolecule Mol){
        CATS2D_7_DL = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_7_DL = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_7_DL = MISSING_VALUE;
            return;
        }
        try {
            ConnAugMatrixCATS = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_7_DL = MISSING_VALUE;
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
                if(i == 7)
                    CATS2D_7_DL = desc[i];
        }
    }

    private void CalculateTopologicalCharge(InsilicoMolecule Mol){
        JGI3 = 0;

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            JGI3 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            JGI3 = MISSING_VALUE;
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

            int[] PathCount = new int[8];
            double[] GGIval = new double[8];
            for (int i=0; i<8; i++) {
                PathCount[i] = 0;
                GGIval[i] = 0;
            }

            for (int i=0; i<nSK; i++)
                for (int j=i; j<nSK; j++) {

                    int CurPath = TopoMat[i][j];
                    if  ((CurPath>0) && (CurPath<=8)) {
                        PathCount[CurPath-1]++;
                        GGIval[CurPath-1] += Math.abs(CTMatrix[i][j]);
                    }
                }

            // Sets descriptors
            for (int i=0; i<8; i++) {
                if(i==2){
                    if (GGIval[i]>0) {
                        JGI3 = GGIval[i] / (double)PathCount[i];
                    } else {
                        JGI3 = 0;
                    }
                }
            }
        } catch (Throwable e) {
            JGI3 = MISSING_VALUE;
        }

    }

    private void CalculateACF(InsilicoMolecule Mol){
        C_026 = 0;
        DescriptorBlock block = new ACF();
        block.Calculate(Mol);
        try {
            C_026 = block.GetByName("C-026").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
            C_026 = MISSING_VALUE;
        }
    }

    private void CalculateFG(InsilicoMolecule Mol){
        nArNH2 = 0;
        DescriptorBlock block = new FunctionalGroups();
        block.Calculate(Mol);
        try {
            nArNH2 = block.GetByName("nArNH2").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
            nArNH2 = MISSING_VALUE;
        }
    }

    private void CalculateConstitutional(InsilicoMolecule Mol){
        Mw = 0; OPerc = 0; nS = 0;
        DescriptorBlock block = new Constitutional();
        block.Calculate(Mol);
        try {
            Mw = block.GetByName("MW").getValue();
            OPerc = block.GetByName("OPerc").getValue();
            nS = block.GetByName("nS").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateMLogP(InsilicoMolecule Mol) {
        MLogP = 0;
        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        try {
            MLogP = descriptorMLogP.MLogP;
        } catch (Exception ex){
            log.warn(ex.getMessage());
            MLogP = MISSING_VALUE;
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
                if (ConnAugMatrixCATS[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrixCATS[i][j]==1) {
                            if (ConnAugMatrixCATS[j][j] == 8) {
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1)
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrixCATS[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrixCATS[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrixCATS[i][j] == 2) && (ConnAugMatrixCATS[j][j] == 8) )
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] > 1.5) ) {
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] != 1) ) {
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
    private double GetResonanceIntegral(IBond bnd) {

        IAtom atA = bnd.getAtom(0);
        IAtom atB = bnd.getAtom(1);
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
