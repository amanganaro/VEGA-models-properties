package insilico.fish_nic.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.fish_nic.descriptors.weights.*;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    private double nCIR = MISSING_VALUE;
    private double nTB = MISSING_VALUE;
    private double nP = MISSING_VALUE;
    private double nCL = MISSING_VALUE;
    private double nR10 = MISSING_VALUE;
    private double TI1 = MISSING_VALUE;
    private double S2K = MISSING_VALUE;
    private double T_N_P = MISSING_VALUE;
    private double T_P_CL = MISSING_VALUE;
    private double T_CL_CL = MISSING_VALUE;
    private double piPC09 = MISSING_VALUE;
    private double PCR = MISSING_VALUE;
    private double Xindex = MISSING_VALUE;
    private double MATS1e = MISSING_VALUE;
    private double GATS7m = MISSING_VALUE;
    private double EEig14ed = MISSING_VALUE;
    private double EEig14dm = MISSING_VALUE;
    private double ESpm1dm = MISSING_VALUE;
    private double GGI4 = MISSING_VALUE;
    private double Seigv = MISSING_VALUE;
    private double MPC2 = MISSING_VALUE;

    private double Mw = 0;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    public double[] GetDescriptors(){
        return new double[]{
                nCIR, nTB, nP, nCL, nR10, TI1, S2K, T_N_P, T_P_CL, T_CL_CL, piPC09, PCR, Xindex, MATS1e, GATS7m, EEig14ed, EEig14dm, ESpm1dm, GGI4, Seigv, MPC2
        };
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateConstitutional(Mol);
        CalculateRings(Mol);
        CalculateTopologicaDistances(Mol);
        CalculateWalkAndPath(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateEdgeAdjacency(Mol);
        CalculateTopologicalCharge(Mol);
        CalculateEigenvalueBased(Mol);
    }

    private void CalculateEigenvalueBased(InsilicoMolecule Mol){

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setSeigv(MISSING_VALUE);
            return;
        }

        double[][] ConnMatrix;
        try {
            ConnMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setSeigv(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        double[] w = VanDerWaals.getWeights(m);
        double refW = VanDerWaals.GetVdWVolume("C");

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

        double SEig=0;
        for (int i=0; i<eigenvalues.length; i++) {
            SEig += eigenvalues[i];
        }

        this.setSeigv(SEig);


    }

    private void CalculateTopologicalCharge(InsilicoMolecule Mol){

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setGGI4(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setGGI4(MISSING_VALUE);
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

            int[] PathCount = new int[5];
            double[] GGIval = new double[5];
            for (int i=0; i<5; i++) {
                PathCount[i] = 0;
                GGIval[i] = 0;
            }

            for (int i=0; i<nSK; i++)
                for (int j=i; j<nSK; j++) {

                    int CurPath = TopoMat[i][j];
                    if  ((CurPath>0) && (CurPath<=5)) {
                        PathCount[CurPath-1]++;
                        GGIval[CurPath-1] += Math.abs(CTMatrix[i][j]);
                    }
                }

            // Sets descriptors
            for (int i=0; i<5; i++) {
                if(i == 3){
                    if (GGIval[i]>0) {
                        this.setGGI4(GGIval[i]);
                    } else {
                       this.setGGI4(0);
                    }
                }

            }

        } catch (Throwable e) {
            this.setGGI4(MISSING_VALUE);
        }


    }

    private void CalculateESPM(InsilicoMolecule Mol){
        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setESpm1dm(MISSING_VALUE);
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            this.setESpm1dm(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setESpm1dm(MISSING_VALUE);
            return;
        }

        Matrix DataMatrix = null;
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

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        // Spectral moment
        for (int i=1; i<=15; i++) {
            double curSM = 0;
            for (int k=(eigenvalues.length-1); k>=0; k--) {
                curSM += Math.pow(eigenvalues[k], (i));
            }
            curSM = Math.log(1 + curSM);

            if(i == 1)
                this.setESpm1dm(curSM);

        }
    }

    private void CalculateEdgeAdjacency(InsilicoMolecule Mol){

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        Matrix DataMatrix;
        // ed
        double[][] EdgeDegreeMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                EdgeDegreeMat[i][j] = EdgeAdjMat[i][j][1];

        DataMatrix = new Matrix(EdgeDegreeMat);

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        // EEig
        for (int i=1; i<=15; i++) {
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                if (i == 14)
                    this.setEEig14ed(eigenvalues[idx]);

        }

        // dm
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
        for (int i=1; i<=15; i++) {
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                if (i == 14)
                    this.setEEig14dm(eigenvalues[idx]);

        }

        // Spectral moment
        for (int i=1; i<=15; i++) {
            double curSM = 0;
            for (int k=(eigenvalues.length-1); k>=0; k--) {
                curSM += Math.pow(eigenvalues[k], (i));
            }
            curSM = Math.log(1 + curSM);

            if(i == 1)
                this.setESpm1dm(curSM);
        }
    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMATS1e(MISSING_VALUE); this.setGATS7m(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setMATS1e(MISSING_VALUE); this.setGATS7m(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = Electronegativity.getWeights(m);
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 1) {
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

        this.setMATS1e(MoranAC);

        w = Mass.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        double GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 7) {

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

        this.setGATS7m(GearyAC);
    }

    private void CalculateWalkAndPath(InsilicoMolecule Mol) {

        List<Integer> PathList = new ArrayList<>();
        for(int i = 1; i<=10; i++)
            PathList.add(i);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            return;
        }

        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();
        double piID=0, TPC=0;
        for (Integer curPath : PathList) {


            double CurPath=0, CurWalk=0, CurPW=0, CurMPC=0;
            int[] AtomWalk = GetAtomsWalks(curPath, AdjMatDbl);
            double[][] AtomPath = GetAtomsPaths(curPath, m);

            for (int i=0; i<nSK; i++) {
                CurWalk += AtomWalk[i];
                CurPath += AtomPath[i][0];
                CurMPC += AtomPath[i][1];
                CurPW += (double)AtomPath[i][0] / (double)AtomWalk[i];
            }

            if (curPath == 1) {
                CurWalk /= 2;
            } else {
                CurWalk = Math.log(1 + CurWalk);
            }
            CurPath /= 2;
            TPC += CurPath;
            CurMPC /= 2.0;
            piID += CurMPC;
            CurMPC = Math.log(1+CurMPC);
            CurPW /= nSK;


            if(curPath == 2)
                this.setMPC2(CurPath);

            if(curPath == 9)
                this.setPiPC09(CurMPC);

        }

        piID = Math.log(1+piID+nSK);
        TPC = Math.log(1+TPC+nSK);

        double PCR = piID / TPC;
        this.setPCR(PCR);

    }

    private void CalculateTopologicaDistances(InsilicoMolecule Mol){
        DescriptorBlock block = new TopologicalDistances();
        block.Calculate(Mol);
        try {
            this.setT_N_P(block.GetByName("T(N..P)").getValue());
            this.setT_CL_CL(block.GetByName("T(Cl..Cl)").getValue());
            this.setT_P_CL(block.GetByName("T(P..Cl)").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setT_N_P(MISSING_VALUE); this.setT_CL_CL(MISSING_VALUE); this.setT_P_CL(MISSING_VALUE);
        }

    }

    private void CalculateRings(InsilicoMolecule Mol){

        try {

            int nSizes = 9;
            int[] RingCount = new int[nSizes];
            int[] RingSize = new int[nSizes];
            for (int i=0; i<nSizes; i++) {
                RingSize[i] = 3 + i;
                RingCount[i] = 0;
            }

            IRingSet allRings = Mol.GetAllRings();
            for (IAtomContainer iAtomContainer : allRings.atomContainers()) {
                IRing ring = (IRing) iAtomContainer;
                for (int i = 0; i < nSizes; i++) {
                    if (ring.getAtomCount() == RingSize[i])
                        RingCount[i]++;
                }
            }

            for (int i=0; i<nSizes; i++)
                if(RingSize[i] == 10)
                    this.setNR10(RingCount[i]);

        } catch (Throwable e) {
            this.setNR10(MISSING_VALUE);
        }
    }

    private void CalculateConstitutional(InsilicoMolecule Mol){

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMw(MISSING_VALUE); this.setNTB(MISSING_VALUE);
            this.setNP(MISSING_VALUE); this.setNCL(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nTotH=0;
            int nC=0, nN=0, nO=0, nP=0, nS=0;
            int nI=0, nF=0, nCl=0, nBr=0, nB=0;
            int nHet=0;
            double mw=0, amw=0, sv=0, mv=0, sp=0, mp=0, se=0, me=0;


            //// Counts on atoms

            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];


                if (CurAt.getSymbol().equalsIgnoreCase("C"))
                    nC++;
                else
                    nHet++;

                if (CurAt.getSymbol().equalsIgnoreCase("P"))
                    nP++;
                if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                    nCl++;
            }

            this.setNP(nP);
            this.setNCL(nCl);
            //// Counts on bonds

            int nArBonds=0, nDblBonds=0, nTrpBonds=0, nMulBonds=0;
            double scbo=0;

            for (int i=0; i<nBO; i++) {

                IBond CurBo = curMol.getBond(i);

                if (CurBo.getFlag(CDKConstants.ISAROMATIC)) {
                    nArBonds++;
                    nMulBonds++;
                    scbo += 1.5;
                } else {
                    if (CurBo.getOrder() == IBond.Order.SINGLE) {
                        scbo++;
                    } else {
                        nMulBonds++;
                        if (CurBo.getOrder() == IBond.Order.DOUBLE) {
                            nDblBonds++;
                            scbo += 2;
                        }
                        if (CurBo.getOrder() == IBond.Order.TRIPLE) {
                            nTrpBonds++;
                            scbo += 3;
                        }
                    }
                }

            }

            this.setNTB(nTrpBonds);

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");
            double[] wVdW = VanDerWaals.getWeights(curMol);
            double HVdW = VanDerWaals.GetVdWVolume("H");
            double[] wPol = Polarizability.getWeights(curMol);
            double HPol = Polarizability.GetPolarizability("H");
            double[] wEl = Electronegativity.getWeights(curMol);
            double HEl = Electronegativity.GetElectronegativity("H");


            for (int i=0; i<nSK; i++) {
                if (wMass[i] == -999)
                    mw = -999;
                if (wVdW[i] == -999)
                    sv = -999;
                if (wPol[i] == -999)
                    sp = -999;
                if (wEl[i] == -999)
                    se = -999;
            }

            for (int i=0; i<nSK; i++) {
                if (mw != -999) {
                    mw += wMass[i];
                    if (H[i]>0) {
                        mw += HMass * H[i];
                    }
                }
                if (sv != -999) {
                    sv += wVdW[i];
                    if (H[i]>0)
                        sv += HVdW * H[i];
                }
                if (sp != -999) {
                    sp += wPol[i];
                    if (H[i]>0)
                        sp += HPol * H[i];
                }
                if (se != -999) {
                    se += wEl[i];
                    if (H[i]>0)
                        se += HEl * H[i];
                }
            }

            if (mw != -999)
                amw = mw/(nSK + nTotH);
            if (sv != -999)
                mv = sv/(nSK + nTotH);
            if (sp != -999)
                mp = sp/(nSK + nTotH);
            if (se != -999)
                me = se/(nSK + nTotH);

            this.setMw(mw);

        } catch (Throwable e) {
            this.setMw(MISSING_VALUE); this.setNTB(MISSING_VALUE);
            this.setNP(MISSING_VALUE); this.setNCL(MISSING_VALUE);
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
