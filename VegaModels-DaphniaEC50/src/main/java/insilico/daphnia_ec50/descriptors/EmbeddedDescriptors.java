package insilico.daphnia_ec50.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.daphnia_ec50.descriptors.weights.*;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;

@Slf4j
public class EmbeddedDescriptors {

    private final int MISSING_VALUE = -999;

    public double EEig15bo = MISSING_VALUE;
    public double ALogP = MISSING_VALUE;
    public double B2_C_O = MISSING_VALUE;
    public double S_106 = MISSING_VALUE;
    public double EEig8dm = MISSING_VALUE;
    public double F4_Cl_Cl= MISSING_VALUE;
    public double GATS1m = MISSING_VALUE;
    public double MATS4p = MISSING_VALUE;
    public double B10_C_N = MISSING_VALUE;
    public double MATS5e = MISSING_VALUE;
    public double Me = MISSING_VALUE;
    public double F10_O_O = MISSING_VALUE;
    public double Mw = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateMe(Mol);
        CalculateACF(Mol);
        CalculateAutocorrelation(Mol);
        CalculateFG(Mol);
        CalculateEAC(Mol);
        CalculateALogP(Mol);
        CalculateEeig15bo(Mol);
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        ALogP = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ALogP = MISSING_VALUE;
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF = new ACF();
        CurACF.Calculate(Mol);

        double[] Frags = CurACF.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            ALogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }


    }

    private void CalculateEeig15bo(InsilicoMolecule Mol){
        EEig15bo = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
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
            if(i==15) {
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    EEig15bo = eigenvalues[idx];
                else
                EEig15bo = 0;
            }

        }


    }

    private void CalculateEAC(InsilicoMolecule Mol) {

        EEig8dm = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            EEig8dm = MISSING_VALUE;
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            EEig8dm = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            EEig8dm = MISSING_VALUE;
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
            int idx = (eigenvalues.length - 1) - (i-1);
            if (idx>=0)
                if(i==8)
                    EEig8dm = eigenvalues[idx];
        }
    }

    private void CalculateFG(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorBlock fg = new TD();
        fg.Calculate(Mol);
        B2_C_O = fg.GetByName("B2(C..O)").getValue();
        F4_Cl_Cl = fg.GetByName("F4(Cl..Cl)").getValue();
        B10_C_N = fg.GetByName("B10(C..N)").getValue();
        F10_O_O = fg.GetByName("F10(O..O)").getValue();
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol){
        GATS1m = 0;
        MATS4p = 0;
        MATS5e = 0;

        int lag_m = 1;
        int lag_p = 4;
        int lag_e = 5;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            GATS1m = MISSING_VALUE;
            MATS4p = MISSING_VALUE;
            MATS5e = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            GATS1m = MISSING_VALUE;
            MATS4p = MISSING_VALUE;
            MATS5e = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w_m = Mass.getWeights(m);
        double[] w_p = Polarizability.getWeights(m);
        double[] w_e = Electronegativity.getWeights(m);

        double wA_m = 0;
        double wA_p = 0;
        double wA_e = 0;
        for (int i=0; i<nSK; i++) {
            wA_m += w_m[i];
            wA_p += w_p[i];
            wA_e += w_e[i];
        }

        wA_m = wA_m / ((double) nSK);
        wA_p = wA_p / ((double) nSK);
        wA_e = wA_e / ((double) nSK);

        double GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_m[i] - wA_m), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_m) {
                    GearyAC += Math.pow((w_m[i] - w_m[j]), 2);
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
        GATS1m = GearyAC;

        double MoranAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_e[i] - wA_e), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_e) {
                    MoranAC += (w_e[i] - wA_e) * (w_e[j] - wA_e);
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


        MATS5e = MoranAC;

        MoranAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_p[i] - wA_p), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_p) {
                    MoranAC += (w_p[i] - wA_p) * (w_p[j] - wA_p);
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


        MATS4p = MoranAC;



    }

    private void CalculateACF(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);
        S_106 = acf.GetByName("S-106").getValue();
    }

    private void CalculateMe(InsilicoMolecule Mol){
        Me = 0; Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            Me = MISSING_VALUE;
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
                } catch (Exception e) {
                    Me = MISSING_VALUE;
                    log.warn(e.getMessage());
                }
                nTotH += H[i];


                if (CurAt.getSymbol().equalsIgnoreCase("C"))
                    nC++;
                else
                    nHet++;

                if (CurAt.getSymbol().equalsIgnoreCase("N"))
                    nN++;
                if (CurAt.getSymbol().equalsIgnoreCase("O"))
                    nO++;
                if (CurAt.getSymbol().equalsIgnoreCase("P"))
                    nP++;
                if (CurAt.getSymbol().equalsIgnoreCase("S"))
                    nS++;
                if (CurAt.getSymbol().equalsIgnoreCase("F"))
                    nF++;
                if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                    nCl++;
                if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                    nBr++;
                if (CurAt.getSymbol().equalsIgnoreCase("I"))
                    nI++;
                if (CurAt.getSymbol().equalsIgnoreCase("B"))
                    nB++;

            }


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
                        Mw += HMass * H[i];
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
                Me = se/(nSK + nTotH);



        } catch (Throwable e) {
            log.warn(e.getMessage());
            Me = MISSING_VALUE;
        }
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
