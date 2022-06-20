package insilico.moa_epa.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.moa_epa.descriptors.weights.*;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    private double MDE_C_22 = MISSING_VALUE;
    private double MATS1m = MISSING_VALUE;
    private double MATS2e = MISSING_VALUE;
    private double GATS2e = MISSING_VALUE;
    private double MDE_C_13 = MISSING_VALUE;
    private double MDE_C_23 = MISSING_VALUE;
    private double BEL1e = MISSING_VALUE;
    private double nR10 = MISSING_VALUE;
    private double GATS6v = MISSING_VALUE;
    private double Gmin = MISSING_VALUE;
    private double SdssNp = MISSING_VALUE;
    private double MDE_C_33 = MISSING_VALUE;
    private double SRW7 = MISSING_VALUE;
    private double SdssC = MISSING_VALUE;
    private double StN = MISSING_VALUE;
    private double MDE_C_34 = MISSING_VALUE;
    private double nX = MISSING_VALUE;
    private double nR5 = MISSING_VALUE;
    private double nR6 = MISSING_VALUE;
    private double StsC = MISSING_VALUE;
    private double SsOm = MISSING_VALUE;
    private double BEH4p = MISSING_VALUE;
    private double Hmin = MISSING_VALUE;
    private double Qsv = MISSING_VALUE;
    private double IVDEM = MISSING_VALUE;
    private double Ms = MISSING_VALUE;
    private double ATS4m = MISSING_VALUE;
    private double GATS8m = MISSING_VALUE;


    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescrpitors(Mol);
    }

    private void CalculateAllDescrpitors(InsilicoMolecule Mol){
        CalculateDistanceEdge(Mol);
        CalculateAutCorrelation(Mol);
        CalculateBEHFilled(Mol);
        CalculateRings(Mol);
        CalculateEStates(Mol);
        CalculateWAP(Mol);
        CalculateInformationContent(Mol);
        CalculateNx(Mol);
    }

    private void CalculateNx(InsilicoMolecule Mol){
        this.setNX(0);
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setNX(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nI=0, nF=0, nCl=0, nBr=0;


            //// Counts on atoms

            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }



                if (CurAt.getSymbol().equalsIgnoreCase("F"))
                    nF++;
                if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                    nCl++;
                if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                    nBr++;
                if (CurAt.getSymbol().equalsIgnoreCase("I"))
                    nI++;


            }


            this.setNX(nI + nF + nCl + nBr);



        } catch (Throwable e) {
            this.setNX(MISSING_VALUE);
        }
    }

    private void CalculateInformationContent(InsilicoMolecule Mol){
        setIVDEM(0);

        int MaxPath = 1;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setIVDEM(MISSING_VALUE);
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
            setIVDEM(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);


        // Information content

        int[] TopoDistFreq = new int[nSK];  // frequencies of topological distances
        double TopoDistFreqSum = 0;
        for (int i=0; i<nSK; i++)
            TopoDistFreq[i] = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i+1; j<nSK; j++) {
                TopoDistFreq[TopoDistMat[i][j]]++;
                TopoDistFreqSum += TopoDistMat[i][j];
            }


        int[] VerDegCount = new int[10];
        for (int i=0; i<10; i++) VerDegCount[i] = 0;
        for (int i=0; i<nSK; i++) {
            VerDegCount[VertexDeg[i]]++;
        }


        double IVDEM = 0;

        for (int i=0; i<10; i++) {
            if (VerDegCount[i]>0)
                IVDEM = IVDEM - ( ( (double)VerDegCount[i] / (double)nSK ) * Log(2, ( (double)VerDegCount[i] / (double)nSK ) ) );
        }
        this.setIVDEM(IVDEM);

    }

    private void CalculateWAP(InsilicoMolecule Mol){
        setSRW7(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setSRW7(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setSRW7(MISSING_VALUE);
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();

        int SRWMaxPath = 8;
        int[] MolSRW = GetSRWToLag(SRWMaxPath, AdjMatDbl);
        for (int i=1; i<=SRWMaxPath; i++)
            if(i==7)
                this.setSRW7(MolSRW[i]);
    }

    private void CalculateEStates(InsilicoMolecule Mol) {
         setGmin(0); setSdssNp(0); setSdssC(0); setStN(0);
         setStsC(0); setSsOm(0); setHmin(0); setQsv(0); setMs(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setGmin(MISSING_VALUE); setSdssNp(MISSING_VALUE); setSdssC(MISSING_VALUE); setStN(MISSING_VALUE);
            setStsC(MISSING_VALUE); setSsOm(MISSING_VALUE); setHmin(MISSING_VALUE); setQsv(MISSING_VALUE); setMs(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            setGmin(MISSING_VALUE); setSdssNp(MISSING_VALUE); setSdssC(MISSING_VALUE); setStN(MISSING_VALUE);
            setStsC(MISSING_VALUE); setSsOm(MISSING_VALUE); setHmin(MISSING_VALUE); setQsv(MISSING_VALUE); setMs(MISSING_VALUE);            return;
        }

        // Calculation
        double SsCl=0, SsF=0, SsBr=0, SsI=0;
        double StN=0, SsssN=0, SsSH=0;
        double SdssC=0, SHssNH=0, SdssNp=0;
        double StsC=0, SsOm=0;

        double Hmax= MISSING_VALUE, Hmin= MISSING_VALUE;
        double Gmax= MISSING_VALUE, Gmin= MISSING_VALUE;

        double Ss = 0, Ms = 0;

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("unable to get H count");
            }

            // formal charge
            int Charge;
            try {
                Charge = curAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }

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

            // Sum of e-states
            Ss += es.getEState()[at];

            // Maximum and minimum Estate/HEstate
            Gmax = (Gmax== Descriptor.MISSING_VALUE) ? es.getEState()[at] : (Math.max(es.getEState()[at], Gmax));
            Gmin = (Gmin== Descriptor.MISSING_VALUE) ? es.getEState()[at] : (Math.min(es.getEState()[at], Gmin));
            if (nH>0) {
                Hmax = (Hmax== Descriptor.MISSING_VALUE) ? es.getHEState()[at] : (Math.max(es.getHEState()[at], Hmax));
                Hmin = (Hmin== Descriptor.MISSING_VALUE) ? es.getHEState()[at] : (Math.min(es.getHEState()[at], Hmin));
            }


            // Halo atoms
            if (curAt.getSymbol().equalsIgnoreCase("Cl")) {
                SsCl += es.getEState()[at];
            }

            if (curAt.getSymbol().equalsIgnoreCase("Br")) {
                SsBr += es.getEState()[at];
            }

            if (curAt.getSymbol().equalsIgnoreCase("F")) {
                SsF += es.getEState()[at];
            }

            if (curAt.getSymbol().equalsIgnoreCase("I")) {
                SsI += es.getEState()[at];
            }


            // C Groups
            if (curAt.getSymbol().equalsIgnoreCase("C")) {

                if ((nBnd == 3) && (nDbl == 1) && (nSng == 2))
                    SdssC += es.getEState()[at];

                if ((nBnd == 2) && (nTri == 1) && (nSng == 1))
                    StsC += es.getEState()[at];

            }


            // C Groups
            if (curAt.getSymbol().equalsIgnoreCase("O")) {

                if ((nBnd == 1) && (nSng == 1) && (Charge == -1))
                    SsOm += es.getEState()[at];

            }


            // N Groups
            if (curAt.getSymbol().equalsIgnoreCase("N")) {

                if ((nBnd == 1) && (nTri == 1))
                    StN += es.getEState()[at];

                if ((nBnd == 3) && (nSng == 3))
                    SsssN += es.getEState()[at];

                if ((nBnd == 2) && (nSng == 2) && (nH == 1))
                    SHssNH += es.getHEState()[at];

                if ((nBnd == 3) && (nSng == 2) && (nDbl == 1) && (Charge == 1))
                    SdssNp += es.getEState()[at];
            }


            // S groups
            if (curAt.getSymbol().equalsIgnoreCase("S")) {

                if ((nBnd == 1) && (nSng == 1) && (nH == 1))
                    SsSH += es.getEState()[at];
            }

        }


        // Qs, Qv and Qsv from TEST software

        double []IS = es.getIS();
        double sumIalk = 0.0D, sumI = 0.0D, sumImax = 0.0D;
        for (int i = 0; i < nSK; i++) {

            double DV = m.getConnectedBondsCount(i);
            double D = DV;

            sumIalk += (DV + 1.0D) / D;
            sumI += IS[i];
            if (DV == 1.0D) {
                sumImax += 8.0D;
            } else if (DV == 2.0D) {
                sumImax += 3.5D;
            } else if (DV == 3.0D) {
                sumImax += 2.0D;
            } else if (DV == 4.0D) {
                sumImax += 1.25D;
            } else {
                sumImax += IS[i];
            }
        }


        double sumIave = (sumIalk + sumImax) / 2.0D;
        double Qsv = (sumIave * sumIalk / (sumI * sumI));


        Ms = Ss / m.getAtomCount();

        this.setMs(Ms);
        this.setSdssC(SdssC);
        this.setStN(StN);
        this.setSdssNp(SdssNp);
        this.setStsC(StsC);
        this.setSsOm(SsOm);
        this.setGmin(Gmin);
        this.setHmin(Hmin);
        this.setQsv(Qsv);




    }

    private void CalculateRings(InsilicoMolecule Mol){
        setNR10(0); setNR5(0); setNR6(0);

        try {

            int nSizes = 9;
            int[] RingCount = new int[nSizes];
            int[] RingSize = new int[nSizes];
            for (int i=0; i<nSizes; i++) {
                RingSize[i] = 3 + i;
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

            for (int i=0; i<nSizes; i++) {
                if(RingSize[i] == 5)
                    this.setNR5(RingCount[i]);
                if(RingSize[i] == 6)
                    this.setNR6(RingCount[i]);
                if(RingSize[i] == 10)
                    this.setNR10(RingCount[i]);

            }

        } catch (Throwable e) {
            setNR10(MISSING_VALUE); setNR5(MISSING_VALUE); setNR6(MISSING_VALUE);
        }

    }

    private void CalculateBEHFilled(InsilicoMolecule Mol){
        setBEL1e(0); setBEH4p(0);

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            setBEL1e(MISSING_VALUE); setBEH4p(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            setBEL1e(MISSING_VALUE); setBEH4p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        // p
        double[] w = Polarizability.getWeights(m);
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
            if(i==3) {
                if (i>(eigenvalues.length-1)) {
                    this.setBEH4p(0);
                } else {
                    this.setBEH4p(eigenvalues[eigenvalues.length-1-i]);
                }
            }

        }


        // e
        BurdenMat = BurdenMatrix.getMatrix(m);
        w = Electronegativity.getWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valL;
            if(i==0){
                if (i>(eigenvalues.length-1)) {
                    this.setBEL1e(0);
                } else {
                    this.setBEL1e(eigenvalues[i]);
                }
            }
        }


    }

    private void CalculateAutCorrelation(InsilicoMolecule Mol){
        setATS4m(0); setGATS8m(0); setGATS6v(0);
        setMATS1m(0); setMATS2e(0); setGATS2e(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setATS4m(MISSING_VALUE); setGATS8m(MISSING_VALUE); setGATS6v(MISSING_VALUE);
            setMATS1m(MISSING_VALUE); setMATS2e(MISSING_VALUE); setGATS2e(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setATS4m(MISSING_VALUE); setGATS8m(MISSING_VALUE); setGATS6v(MISSING_VALUE);
            setMATS1m(MISSING_VALUE); setMATS2e(MISSING_VALUE); setGATS2e(MISSING_VALUE);
            return;
        }

        // !!! in origine era usata la topological matrix del cdk

        int nSK = m.getAtomCount();


        //m
        double[] w = Mass.getWeights(m);
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

        this.setMATS1m(MoranAC);


        double AC=0;
        for (int i=0; i<nSK; i++) {
            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 4) {
                    AC += w[i] * w[j];
                }
        }

        AC /= 2.0;
        AC = Math.log(1 + AC);

        this.setATS4m(AC);


        double  GearyAC=0;
        denom = 0; delta = 0;
        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 8) {
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

        this.setGATS8m(GearyAC);

        // v
        w = VanDerWaals.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        GearyAC=0;
        denom = 0; delta = 0;
        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 6) {
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

        this.setGATS6v(GearyAC);

        // e
        w = Electronegativity.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        GearyAC=0;
        denom = 0; delta = 0;
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

        this.setGATS2e(GearyAC);

        MoranAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 2) {
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

        this.setMATS2e(MoranAC);
    }

    private void CalculateDistanceEdge(InsilicoMolecule Mol){
        setMDE_C_13(0); setMDE_C_22(0); setMDE_C_23(0);
        setMDE_C_33(0); setMDE_C_34(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setMDE_C_13(MISSING_VALUE); setMDE_C_22(MISSING_VALUE); setMDE_C_23(MISSING_VALUE);
            setMDE_C_33(MISSING_VALUE); setMDE_C_34(MISSING_VALUE);
            return;
        }

        // Get matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setMDE_C_13(MISSING_VALUE); setMDE_C_22(MISSING_VALUE); setMDE_C_23(MISSING_VALUE);
            setMDE_C_33(MISSING_VALUE); setMDE_C_34(MISSING_VALUE);
            return;
        }

        // Get VD
        int[] VD = VertexDegree.getWeights(m, true);

        int nSK = m.getAtomCount();

//        ArrayList<String> ATList = BuildAtomTypeList();
        ArrayList<Integer> VDList = new ArrayList<>();
        VDList.add(1);
        VDList.add(2);
        VDList.add(3);
        VDList.add(4);

            for (Integer curVD : VDList)
                for (Integer curVD2 : VDList) {

                    double n = 0;
                    double d = 1;

                    for (int i=0; i<nSK; i++)
                        if ("C".equalsIgnoreCase(m.getAtom(i).getSymbol()))
                            for (int j=i+1; j<nSK; j++)
                                if ("C".equalsIgnoreCase(m.getAtom(j).getSymbol()))
                                    if  ( ((VD[i]==curVD) && (VD[j]==curVD2)) || ((VD[i]==curVD2) && (VD[j]==curVD)) ) {
                                        n++;
                                        d *= TopoMatrix[i][j];
                                    }

                    double res = 0;
                    if (n > 0) {
                        d = Math.pow(d,1.0/(2.0*n));
                        res = n/(d*d);
                    }

                    if(curVD == 1 && curVD2 == 3)
                        setMDE_C_13(res);
                    if(curVD == 2)
                    {
                        if(curVD2 == 2)
                            setMDE_C_22(res);
                        if(curVD2 == 3)
                            setMDE_C_23(res);
                    }
                    if (curVD == 3){
                        if (curVD2 == 3)
                            setMDE_C_33(res);
                        if (curVD2 == 4)
                            setMDE_C_34(res);
                    }
                }

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

    private static double Log(int base,double x) {
        double Logbx = Math.log10(x)/Math.log10((double)base);
        return Logbx;
    }

}
