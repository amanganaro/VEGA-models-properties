package insilico.fathead_epa.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.fathead_epa.descriptors.weights.*;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.*;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

@Log4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    private double SsssN = MISSING_VALUE;
    private double MDE_C_33 = MISSING_VALUE;
    private double MDE_O_11 = MISSING_VALUE;
    private double BEH2m = MISSING_VALUE;
    private double nDblBo = MISSING_VALUE;
    private double nS = MISSING_VALUE;
    private double nR9 = MISSING_VALUE;
    private double ATS5p = MISSING_VALUE;
    private double MATS7e = MISSING_VALUE;
    private double GATS3e = MISSING_VALUE;
    private double SRW3 = MISSING_VALUE;
    private double ALogP = MISSING_VALUE;
    private double MW = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateEStates(Mol);
        CalculateDistanceEdge(Mol);
        CalculateBE(Mol);
        CalculateConstitutional(Mol);
        CalculateRings(Mol);
        CalculateWAP(Mol);
        CalculateAutocorrelation(Mol);
        CalculateALogP(Mol);
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        setALogP(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setALogP(MISSING_VALUE);
            return;
        }

        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);

        double LogP = 0;
        double[] Frags = acf.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }

        this.setALogP(LogP);
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol) {
        setATS5p(0); setMATS7e(0); setGATS3e(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setATS5p(MISSING_VALUE); setMATS7e(MISSING_VALUE); setGATS3e(MISSING_VALUE);
            return;
        }

        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setATS5p(MISSING_VALUE); setMATS7e(MISSING_VALUE); setGATS3e(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        //p
        double[] w = Polarizability.getWeights(m);
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double AC=0;

        for (int i=0; i<nSK; i++) {
            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 5) {
                    AC += w[i] * w[j];
                }
        }

        // AC transformed in log form
        AC /= 2.0;
        AC = Math.log(1 + AC);

        this.setATS5p(AC);

        // e
        w = Electronegativity.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 7) {
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

        this.setMATS7e(MoranAC);

        w = Electronegativity.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 3) {
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


        this.setGATS3e(GearyAC);

    }

    private void CalculateWAP(InsilicoMolecule Mol) {
        this.setSRW3(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setSRW3(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setSRW3(MISSING_VALUE);
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();
        double piID=0, TPC=0;

        piID = Math.log(1+piID+nSK);
        TPC = Math.log(1+TPC+nSK);

        double PCR = piID / TPC;

//        int SRWMaxPath = pathsList.get(pathsList.size()-1);
        int[] MolSRW = GetSRWToLag(4, AdjMatDbl);
        for (int i=1; i<=10; i++)
            if(i==3)
                setSRW3(MolSRW[i]);
    }

    private void CalculateRings(InsilicoMolecule Mol){
        this.setNR9(0);

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
                if(RingSize[i] == 9)
                    this.setNR9(RingCount[i]);

        } catch (Throwable e) {
            this.setNR9(MISSING_VALUE);
        }
    }

    private void CalculateConstitutional(InsilicoMolecule Mol) {

        DescriptorBlock block = new Constitutional();
        block.Calculate(Mol);

        try {
            MW = block.GetByName("MW").getValue();
            nDblBo = block.GetByName("nDB").getValue();
            nS = block.GetByName("nS").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
            MW = MISSING_VALUE; nDblBo = MISSING_VALUE; nS = MISSING_VALUE;
        }

    }

    private void CalculateBE(InsilicoMolecule Mol){
        BEH2m = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH2m = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            BEH2m = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w = Mass.getWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
                valL = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
                valL = eigenvalues[i];
            }
            if (i == 1)
                this.setBEH2m(valH);
        }
    }

    private void CalculateDistanceEdge(InsilicoMolecule Mol) {
        MDE_C_33 = 0; MDE_O_11 = 0;
        short vd_c = 3;
        short vd_o = 1;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            MDE_C_33 = MISSING_VALUE; MDE_O_11 = MISSING_VALUE;
            return;
        }

        // Get matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MDE_C_33 = MISSING_VALUE; MDE_O_11 = MISSING_VALUE;
            return;
        }

        // Get VD
        int[] VD = VertexDegree.getWeights(m, true);

        int nSK = m.getAtomCount();

        double n = 0;
        double d = 1;

        for (int i=0; i<nSK; i++)
            if ("C".equalsIgnoreCase(m.getAtom(i).getSymbol()))
                for (int j=i+1; j<nSK; j++)
                    if ("C".equalsIgnoreCase(m.getAtom(j).getSymbol()))
                        if  ( ((VD[i]==vd_c) && (VD[j]==vd_c)) || ((VD[i]==vd_c) && (VD[j]==vd_c)) ) {
                            n++;
                            d *= TopoMatrix[i][j];
                        }

        double res = 0;
        if (n > 0) {
            d = Math.pow(d,1.0/(2.0*n));
            res = n/(d*d);
        }
        this.setMDE_C_33(res);

        n = 0;
        d = 1;

        for (int i=0; i<nSK; i++)
            if ("O".equalsIgnoreCase(m.getAtom(i).getSymbol()))
                for (int j=i+1; j<nSK; j++)
                    if ("O".equalsIgnoreCase(m.getAtom(j).getSymbol()))
                        if  ( ((VD[i]==vd_o) && (VD[j]==vd_o)) || ((VD[i]==vd_o) && (VD[j]==vd_o)) ) {
                            n++;
                            d *= TopoMatrix[i][j];
                        }

        res = 0;
        if (n > 0) {
            d = Math.pow(d,1.0/(2.0*n));
            res = n/(d*d);
        }
        this.setMDE_O_11(res);
    }

    private void CalculateEStates(InsilicoMolecule Mol){

        SsssN = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SsssN = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            SsssN = MISSING_VALUE;
            return;
        }

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count bonds
            int nBnd=0, nSng = 0, nDbl = 0, nTri = 0, nAr=0;
            for (IBond b : m.getConnectedBondsList(curAt)) {
                if (b.getFlag(CDKConstants.ISAROMATIC)) {
                    nAr++;
                    nBnd++;
                    continue;
                }
                if (b.getOrder() == IBond.Order.SINGLE) {
                    nSng++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.DOUBLE) {
                    nDbl++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nTri++;
                    nBnd++;
                }
            }

            // N Groups
            if (curAt.getSymbol().equalsIgnoreCase("N")) {

                if ((nBnd == 3) && (nSng == 3))
                    SsssN += es.getEState()[at];
            }

        }

        this.setSsssN(SsssN);


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

}







