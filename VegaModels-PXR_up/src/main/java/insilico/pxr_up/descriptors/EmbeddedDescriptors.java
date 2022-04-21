package insilico.pxr_up.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.basic.WeightsVanDerWaals;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsQuantumNumber;
import insilico.core.descriptor.blocks.weights.other.WeightsValenceVertexDegree;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.pxr_up.descriptors.weights.*;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private final int MISSING_VALUE = -999;

    public double ChiA_Dz_Z = MISSING_VALUE;
    public double Eta_alpha_A = MISSING_VALUE;
    public double GGI3 = MISSING_VALUE;
    public double MAXDN = MISSING_VALUE;
    public double MLOGP2 = MISSING_VALUE;
    public double P_VSA_ppp_D = MISSING_VALUE;
    public double P_VSA_ppp_L = MISSING_VALUE;
    public double SdssC = MISSING_VALUE;
    public double SpMax_B_m = MISSING_VALUE;
    public double SpMax_B_s = MISSING_VALUE;
    public double SpMax4_Bh_m = MISSING_VALUE;
    public double SpMin2_Bh_m = MISSING_VALUE;
    public double SpPosA_B_v = MISSING_VALUE;
    public double Wap = MISSING_VALUE;
    public double X5v = MISSING_VALUE;

    public double ExperimentalValue = MISSING_VALUE;

    private IAtomContainer CurMol;
    private int nAtoms;
    private double[][] ConnMatrix;
    private int[][] TopoMatrix;
    private RingSet MolRings;

    public double[] getDescriptors(){
        return new double[]{
                ChiA_Dz_Z, Eta_alpha_A, GGI3, MAXDN, MLOGP2, P_VSA_ppp_D, P_VSA_ppp_L, SdssC, SpMax_B_m, SpMax_B_s, SpMax4_Bh_m, SpMin2_Bh_m, SpPosA_B_v, Wap, X5v
        };
    }

    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        if(fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateChiDzZ(mol);
        CalculateEtaAlphaA(mol);
        CalculateGGI(mol);
        CalculateMaxDN(mol);
        CalculateMLOGP2(mol);
        CalculatePVSA(mol);
        CalculateSdssC(mol);
        CalculateSpMatrices(mol);
        CalculateSpBurden(mol);
        CalculateWap(mol);
        CalculateX5v(mol);
    }

    private void CalculateChiDzZ(InsilicoMolecule mol) {
        IAtomContainer curMol = null;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
        }


        try {
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = new double[nSK][nSK];

            double[][][] BarMat = mol.GetMatrixBarysz();
            int BarLayer = 0;
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++)
                    Mat[i][j] = BarMat[i][j][BarLayer];

//            matDesc.Calculate(Mat, nSK, nBo, AdjMatrix);
            // row sums
            double[] VS = new double[nSK];
            for (int i=0; i<nSK; i++) {
                VS[i] = 0;
                for (int j=0; j<nSK; j++)
                    VS[i] += Mat[i][j];
            }

            double Chi = 0;

            for (int i=0; i<(nSK-1); i++)
                for (int j=(i+1); j<nSK; j++) {
                    if (AdjMatrix[i][j] > 0)
                        Chi += Math.pow( (VS[i] * VS[j]) , -0.5);
                }

            ChiA_Dz_Z = Chi / (double)nBo;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }




    }

    private void CalculateX5v(InsilicoMolecule mol) {

        IAtomContainer m = null;
        int MaxPath = 5;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
        }

        // Gets matrices
        double[][] ConnAugMatrix = new double[0][];
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());

        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);
        double[] curDescX = new double[MaxPath+1];
        double[] curDescXv = new double[MaxPath+1];
        double[] curDescXsol = new double[MaxPath+1];

        // checks for missing weights
        for (int i=0; i<Qnumbers.length; i++)
            if (Qnumbers[i] == -999)
                return;
        for (int i=0; i<ValenceVD.length; i++)
            if (ValenceVD[i] == -999)
                return;

        // clears VD matrix from linked F
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ((ConnAugMatrix[i][j]>0) && (ConnAugMatrix[j][j]==9))
                    VD[i]--;
            }


        for (int k=0; k<MaxPath; k++) {
            curDescX[k] = 0;
            curDescXv[k] = 0;
            curDescXsol[k] = 0;
        }

        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account

            // path 0
            curDescX[0] += Math.pow(VD[i], -0.5);
            curDescXv[0] += Math.pow(ValenceVD[i], -0.5);
            curDescXsol[0] += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);

            // path 1 - MaxPath
            for (int path=1; path<(MaxPath+1); path++) {

                if (curDescX[path] == -999) continue;

                IAtom at =  m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (int k=0; k<CurPaths.size(); k++) {
                    double prodX = 1;
                    double prodXv = 1;
                    int prodQuantum = 1;
                    List<IAtom> CurPath = CurPaths.get(k);
                    for (int kk=0; kk<CurPath.size(); kk++) {
                        int atIdx = m.getAtomNumber(CurPath.get(kk));
//                        if (ConnMatrix[atIdx][atIdx] == 9)
//                            continue; // F not taken into account
                        prodX *= VD[atIdx];
                        prodXv *= ValenceVD[atIdx];
                        prodQuantum *= Qnumbers[atIdx];
                    }
                    curDescX[path] += Math.pow(prodX, -0.5);
                    curDescXv[path] += Math.pow(prodXv, -0.5);
                    curDescXsol[path] += (1.00 / Math.pow(2.00, (double)(path + 1))) *
                            ((double)prodQuantum) * Math.pow(prodX, -0.5);
                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath+1); i++) {
            curDescX[i] /= 2;
            curDescXv[i] /= 2;
            curDescXsol[i] /= 2;
        }

        // Sets descriptors
        for (int i=0; i<(MaxPath+1); i++) {
            if(i == 5)
                X5v = curDescXv[i];
        }

    }

    private void CalculateWap(InsilicoMolecule mol) {
        try {

            MoleculePaths paths = new MoleculePaths(mol);

            Wap = 0;
            for (int i=0; i<paths.TotalPC.length; i++)
                Wap += (i+1) * paths.TotalPC[i];

        } catch (GenericFailureException e) {
            log.warn("Error while calculating paths, unable to calculate: " + mol.GetSMILES());
        }
    }

    private void CalculateSpBurden(InsilicoMolecule mol) {
        IAtomContainer m = null;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
        }

        // Gets matrix
        double[][] BurdenMat = new double[0][];
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
        }

        int nSK = m.getAtomCount();
        iBasicWeight weight = new WeightsMass();
        double[] w = weight.getScaledWeights(m);

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
            if (i > (eigenvalues.length - 1)) {
                valH = 0;
                valL = 0;
            } else {
                if (eigenvalues[eigenvalues.length - 1 - i] > 0)
                    valH = eigenvalues[eigenvalues.length - 1 - i];
                else
                    valH = 0;
                if (eigenvalues[i] < 0)
                    valL = Math.abs(eigenvalues[i]);
                else
                    valL = 0;
            }

            if (i == 1)
                SpMin2_Bh_m = valL;

            if (i == 3)
                SpMax4_Bh_m = valH;
        }
    }

    private void CalculateSpMatrices(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure for: " + mol.GetSMILES());
            return;
        }

        try {
            // Adj matrix available for calculations
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = mol.GetMatrixBurden();
            double[] w = new double[nSK];

            // weight m
            w = (new WeightsMass()).getScaledWeights(curMol);
            for (int i=0; i<nSK; i++)
                Mat[i][i] = w[i];

            double[] eigenvalues;
            Matrix DataMatrix = new Matrix(Mat);

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("Invalid structure for: " + e.getMessage());
                return;
            }

            double EigMax = eigenvalues[0];

            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
            }

            SpMax_B_m = EigMax;

            // s
            w = (new WeightsIState()).getWeights(curMol, false);
            Mat = mol.GetMatrixBurden();
            for (int i=0; i<nSK; i++)
                Mat[i][i] = w[i];
            DataMatrix = new Matrix(Mat);

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("Unable to calculate eigenvalues: " + e.getMessage());
                return;
            }

            EigMax = eigenvalues[0];

            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
            }

            SpMax_B_s = EigMax;

            w = (new WeightsVanDerWaals()).getScaledWeights(curMol);
            Mat = mol.GetMatrixBurden();
            for (int i=0; i<nSK; i++)
                Mat[i][i] = w[i];
            DataMatrix = new Matrix(Mat);

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("Unable to calculate eigenvalues: " + e.getMessage());
            }

            double SpPos = 0;

            for (double val : eigenvalues) {
                if (val > 0)
                    SpPos += val;
            }

            SpPosA_B_v = SpPos / (double) nSK;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateSdssC(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(mol.GetStructure());
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        // Calculation
        double SdssC=0;


        double Hmax=Descriptor.MISSING_VALUE, Hmin=Descriptor.MISSING_VALUE;
        double Gmax=Descriptor.MISSING_VALUE, Gmin=Descriptor.MISSING_VALUE;

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
            Gmax = (Gmax==Descriptor.MISSING_VALUE) ? es.getEState()[at] : (Math.max(es.getEState()[at], Gmax));
            Gmin = (Gmin==Descriptor.MISSING_VALUE) ? es.getEState()[at] : (Math.min(es.getEState()[at], Gmin));
            if (nH>0) {
                Hmax = (Hmax==Descriptor.MISSING_VALUE) ? es.getHEState()[at] : (Math.max(es.getHEState()[at], Hmax));
                Hmin = (Hmin==Descriptor.MISSING_VALUE) ? es.getHEState()[at] : (Math.min(es.getHEState()[at], Hmin));
            }


            // Halo atoms


            // C Groups
            if (curAt.getSymbol().equalsIgnoreCase("C")) {

                if ((nBnd == 3) && (nDbl == 1) && (nSng == 2))
                    SdssC += es.getEState()[at];

            }
        }

        this.SdssC = SdssC;
    }

    private void CalculatePVSA(InsilicoMolecule mol) {


        // Calculate PPP before everything else, as the original H-depleted molecule is used
        IAtomContainer curMolNoH = null;
        double[][] ConnAugMatrixNoH = null;
        try {
            curMolNoH = mol.GetStructure();
            ConnAugMatrixNoH = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn("Invalid structure for: " +mol.GetSMILES());
        }

        int nSKnoH = curMolNoH.getAtomCount();

//        Cats2D cats = new Cats2D();
        ArrayList<String>[] AtomTypesOnMolWithoutH = setCatsAtomType(curMolNoH, ConnAugMatrixNoH);

        IAtomContainer m = null;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure for: " + mol.GetSMILES());
        }

        int nSK = m.getAtomCount();
        int nBO = m.getBondCount();

        GhoseCrippenACF GC = new GhoseCrippenACF(m, true);
        int[] ACF = GC.GetACF();

        double[] VSA = new double[nSK];
        double[] PartialVSA = new double[nSK];
        double[] IdealBondLen = new double[nBO];
        double[] VdWRadius = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = m.getAtom(i);

            VSA[i] = 0;
            PartialVSA[i] = 0;

            // R - van der waals radius
            VdWRadius[i] = GetVdWRadius(m, at);
            if (VdWRadius[i] == Descriptor.MISSING_VALUE) {
                log.warn("Missing VSA");
                VSA[i] = Descriptor.MISSING_VALUE;
                PartialVSA[i] = Descriptor.MISSING_VALUE;
            }
        }

        for (int i=0; i<nBO; i++) {
            IBond bo = m.getBond(i);

            // refR - reference bond length
            double refR = this.GetRefBondLength(bo.getAtom(0), bo.getAtom(1));
            if (refR == Descriptor.MISSING_VALUE) {
                log.warn("Missing radius");
            }

            // c - correction
            double c = 0;
            if (refR != 0) {
                double bnd = MoleculeUtilities.Bond2Double(bo);
                if (bnd == 1.5) c = 0.1;
                if (bnd == 2) c = 0.2;
                if (bnd == 3) c = 0.3;
            }

            IdealBondLen[i] = refR - c;
        }

        for (int b=0; b<nBO; b++) {
            IBond bo = m.getBond(b);

            int at1 = m.indexOf(bo.getAtom(0));
            int at2 = m.indexOf(bo.getAtom(1));
            if ( (IdealBondLen[b] == 0) || (VdWRadius[at1] == Descriptor.MISSING_VALUE) || (VdWRadius[at2] == Descriptor.MISSING_VALUE) ) {
                PartialVSA[at1] = Descriptor.MISSING_VALUE;
                PartialVSA[at2] = Descriptor.MISSING_VALUE;
                VSA[at1] = Descriptor.MISSING_VALUE;
                VSA[at2] = Descriptor.MISSING_VALUE;
                continue;
            }

            double curDistance;
            double diffRadius = Math.abs(VdWRadius[at1]-VdWRadius[at2]);
            if (diffRadius > IdealBondLen[b])
                curDistance = diffRadius;
            else
                curDistance = IdealBondLen[b];

            if (curDistance > (VdWRadius[at1]+VdWRadius[at2]))
                curDistance = VdWRadius[at1]+VdWRadius[at2];

            PartialVSA[at1] = PartialVSA[at1] + ( Math.pow(VdWRadius[at2],2) - Math.pow((VdWRadius[at1] - curDistance), 2) ) / curDistance;
            PartialVSA[at2] = PartialVSA[at2] + ( Math.pow(VdWRadius[at1],2) - Math.pow((VdWRadius[at2] - curDistance), 2) ) / curDistance;
        }

        for (int i=0; i<nSK; i++) {
            if (VSA[i] != Descriptor.MISSING_VALUE)
                VSA[i] = 4 * Math.PI * Math.pow(VdWRadius[i], 2) - (Math.PI * VdWRadius[i] * PartialVSA[i]);
            if (VSA[i] < 0)
                VSA[i] = 0;
        }

        int bins = 5;
        for (int b=1; b<=bins; b++) {

            double PVSA = 0;



            // special case for PPP
            // in the new mol the H are all added after the original atoms, so the first nSK atoms
            // and calculated VSA are the same both for the original mol and for the H-filled one

            String curAtomType = PPP_TYPES[b-1];

            for (int i = 0; i < nSKnoH; i++) {

                for (String atType : AtomTypesOnMolWithoutH[i]){
                    if (atType.equalsIgnoreCase(curAtomType)) {
                        PVSA += VSA[i];
                        break;
                    }
                }
            }

            if(curAtomType.equalsIgnoreCase("D")){
                P_VSA_ppp_D = PVSA;
            }

            if(curAtomType.equalsIgnoreCase("L")){
                P_VSA_ppp_L = PVSA;
            }

        }
    }

    private void CalculateMLOGP2(InsilicoMolecule mol) {

        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
        }

        // Gets matrices
        try {
            ConnMatrix = mol.GetMatrixConnectionAugmented();
            TopoMatrix = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
        }

        // !!! in origine usava topo matrix di CDK

        double CX = 0;
        int NO = 0;
        int PRX = 0;
        int UB = 0;
        int HB = 0;
        int POL = 0;
        double AMP = 0;
        int ALK = 0;
        int RNG = 0;
        double QN = 0;
        int NO2 = 0;
        double NCS = 0;
        int BLM = 0;

        double LogP = 0;

        nAtoms = CurMol.getAtomCount();
        try {
            MolRings = mol.GetSSSR();
        } catch (InvalidMoleculeException e) {
            return;
        }


        // Atom Counters

        int nC=0, nN=0, nO=0, nF=0, nP=0, nS=0, nCl=0, nBr=0, nI=0;
        int nNO2=0, nNamm=0, nNoxide=0;
        for (int i=0; i<nAtoms; i++) {
            int atom = (int)ConnMatrix[i][i];
            switch (atom) {
                case 6: nC++; break;
                case 7: nN++; break;
                case 8: nO++; break;
                case 9: nF++; break;
                case 15: nP++; break;
                case 16: nS++; break;
                case 17: nCl++; break;
                case 35: nBr++; break;
                case 53: nI++; break;
            }

            // N centered groups
            if (atom == 7) {

                int nOdbl = 0;
                int nOminus = 0;
                int nVD = 0;
                int nCoxide = 0;
                double totBondOrder = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (j == i) continue;
                    if (ConnMatrix[i][j] > 0) {
                        nVD++;
                        totBondOrder += ConnMatrix[i][j];
                        if  (ConnMatrix[j][j] == 6)
                            nCoxide++;
                        if ((ConnMatrix[i][j] == 2) && (ConnMatrix[j][j] == 8))
                            nOdbl++;
                        if ((ConnMatrix[i][j] == 1) && (ConnMatrix[j][j] == 8)) {
                            int OVD = 0;
                            for (int k = 0; k < nAtoms; k++) {
                                if (k == j) continue;
                                if (ConnMatrix[j][k] > 0) OVD++;
                            }
                            if (OVD == 1)
                                nOminus++;
                        }
                    }
                }
                int nH;
                try {
                    nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }

                // count quaternary N (or N oxide with VD = 4)
                if ((nVD+nH) == 4) {
                    if ((nOdbl == 1) && (nCoxide + nH == 3))
                        nNoxide++;
                    else
                        nNamm++;
                }

                // count N oxide with VD = 2
                if ( (nVD == 2) && (totBondOrder == 4) )
                    if (nOdbl > 0)
                        nNoxide++;

                // N oxide in pyridine type structure
                if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC))
                    if ( (nVD == 3) && (nOdbl == 1) && (nCoxide == 2) )
                        nNoxide++;

                // in Dragon seems to ignore other possible N oxides ? possible error

                // count NO2 groups
                if (InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION) {
                    // due to preprocessing, nitro groups are  set in O=N=O form
                    if ((nVD == 3) && (nOdbl == 2) )
                        nNO2++;
                } else {
                    // due to preprocessing, nitro groups are  set in O=[N+][O-] form
                    if ((nVD == 3) && (nOdbl == 1) && (nOminus == 1))
                        nNO2++;
                }
            }
        }


        //// 1 - CX

        CX = nC + 0.5 * nF + nCl + 1.5 * nBr + 2 * nI;


        //// 2 - NO

        NO = nN + nO;


        //// 3 - PRX

        PRX = Calculate_PRX();


        //// 4 - UB

        for (IBond bnd : CurMol.bonds()) {
            // new CDK correctly parses aromatic bonds to kekule structure
            // indeed their order is set correctly to single or double
            if ( (bnd.getOrder().numeric() == 2) || (bnd.getOrder().numeric() == 3) )
                UB++;
        }

        // don't count double bonds in NO2 (consider one dbl for each NO2 group with standard normalization)
        int DblBondsInNO2 = 1;
        if (InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION)
            DblBondsInNO2 = 2;
        UB -= nNO2 * DblBondsInNO2;

        // Correction for aromatic N like in pyridine oxide
        for (int i=0; i<nAtoms; i++) {
            if ( (ConnMatrix[i][i] == 7) && (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) ) {
                // TODO: ma cosa fa in dragon??
            }
        }


        //// 5. HB
        HB = Check_HB() ? 1 : 0;


        //// 6. POL
        POL = Calculate_POL();
        if (POL > 4)
            POL = 4;


        //// 7. AMP
        AMP = Calculate_AMP(CurMol);


        //// 8. ALK
        ALK = Calculate_ALK(CurMol);


        //// 9. RNG
        RNG = Calculate_RNG(CurMol);


        //// 10. QN
        QN = nNamm + 0.5 * nNoxide;


        //// 11. NO2
        NO2 = nNO2;


        //// 12. NCS
        NCS = Check_NCS();


        //// 13. BLM

        for (int i=0; i<nAtoms; i++) {

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("N")) {

                if (BLM == 0)
                    BLM = Check_BLM(CurMol, i);

            }

        }



        // Calculation of LogP
        LogP = -1.041 + 1.244 * Math.pow(CX,0.6) - 1.017 * Math.pow(NO, 0.9) +
                0.406 * PRX + -0.145 * Math.pow(UB, 0.8) + 0.511 * HB +
                0.268 * POL - 2.215 * AMP + 0.912 * ALK -0.392 * RNG -
                3.684 * QN + 0.474 * NO2 + 1.582 * NCS + 0.773 * BLM;

        MLOGP2 = Math.pow(LogP, 2);
    }

    private void CalculateMaxDN(InsilicoMolecule mol) {

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        // Get weights
        EState ES;
        try {
            ES = new EState(mol.GetStructure());
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }
        double[] IStates = ES.getIS();

        // Get matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }


        try {

            int nSK = curMol.getAtomCount();
            double maxDN = 0, maxDP = 0;

            for (int i=0; i<nSK; i++) {
                double Delta = 0;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if ( (IStates[i]!= Descriptor.MISSING_VALUE) && (IStates[j]!= Descriptor.MISSING_VALUE) )
                        Delta += (IStates[i]-IStates[j]) / Math.pow(TopoMat[i][j]+1, 2);
                }
                if (Delta<0) {
                    if (Delta<maxDN)
                        maxDN = Delta;
                } else {
                    if (Delta>maxDP)
                        maxDP = Delta;
                }
            }


            MAXDN = Math.abs(maxDN);

        } catch (Throwable e) {
            log.warn(e.getMessage());
        }

    }

    private void CalculateGGI(InsilicoMolecule mol) {
        int MAXPATH = 8;

        IAtomContainer curMol = null;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure for: " + mol.GetSMILES());
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
        }

        try {

            int nSK = curMol.getAtomCount();
            int[][] AdjMat = AdjacencyMatrix.getMatrix(curMol);
            int[] VertexDegrees = WeightsVertexDegree.getWeights(curMol, false);

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

            int[] PathCount = new int[MAXPATH];
            double[] GGIval = new double[MAXPATH];
            for (int i=0; i<MAXPATH; i++) {
                PathCount[i] = 0;
                GGIval[i] = 0;
            }

            for (int i=0; i<nSK; i++)
                for (int j=i; j<nSK; j++) {

                    int CurPath = TopoMat[i][j];
                    if  ((CurPath>0) && (CurPath<=MAXPATH)) {
                        PathCount[CurPath-1]++;
                        GGIval[CurPath-1] += Math.abs(CTMatrix[i][j]);
                    }
                }

            // Sets descriptors
            for (int i=0; i<MAXPATH; i++) {
                if (GGIval[i]>=0) {
                    if(i == 2)
                        GGI3 = GGIval[i];
                }
            }


        } catch (Throwable e) {
            log.warn(e.getMessage());
        }


    }

    private void CalculateEtaAlphaA(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure, unable to calculate: " + mol.GetSMILES());
            return;
        }

        // Gets matrices
        double[][] AdjConnMatrix;
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(curMol);
            AdjConnMatrix = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        RingSet MolRings;
        try {
            MolRings = mol.GetSSSR();
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = curMol.getAtomCount();

        //// Initial parameters

        int[] AtomicNumber = new int[nSK];
        int[] ElectronNumber = new int[nSK]; // valence electron number
        int[] QuantumNumber = new int[nSK];
        int[] VertexDegree = WeightsVertexDegree.getWeights(curMol, true); // include H in VD
        boolean[] isAromatic = new boolean[nSK];

        WeightsValenceVertexDegree wVVD = new WeightsValenceVertexDegree();
        WeightsQuantumNumber wQN = new WeightsQuantumNumber();
        for (int i=0; i<nSK; i++) {
            AtomicNumber[i] = curMol.getAtom(i).getAtomicNumber();
            ElectronNumber[i] = wVVD.GetValenceElectronsNumber(curMol.getAtom(i).getSymbol());
            QuantumNumber[i] = wQN.getWeight(curMol.getAtom(i).getSymbol());
            if ( (ElectronNumber[i] == Descriptor.MISSING_VALUE) || (QuantumNumber[i] == Descriptor.MISSING_VALUE) ) {
                log.warn("weight missing for atom n. " + (i+1) + ": " + curMol.getAtom(i).getSymbol());
                return;
            }
            isAromatic[i] = curMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);
        }

        //// Parameters for ETA indices
        double[] alpha = new double[nSK];

        for (int i=0; i<nSK; i++) {
            alpha[i] = 0;
        }

        for (int i=0; i<nSK; i++) {
            alpha[i] = (double)(AtomicNumber[i] - ElectronNumber[i]) / ElectronNumber[i] * (1.0 / (QuantumNumber[i] - 1));
        }



        double eta_alpha = 0;

        for (int i=0; i<nSK; i++) {
            eta_alpha += alpha[i];
        }

        Eta_alpha_A = eta_alpha / (double) nSK;
    }

    private void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {
        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-PXR_up/src/main/resources/data/PXR_up.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            br.readLine();
            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[2].trim()).GetSMILES())) {

                    ExperimentalValue = Double.parseDouble(lineArray[8]);

                    ChiA_Dz_Z = Double.parseDouble(lineArray[15]);
                    Eta_alpha_A = Double.parseDouble(lineArray[16]);
                    GGI3 = Double.parseDouble(lineArray[17]);
                    MAXDN = Double.parseDouble(lineArray[18]);
                    MLOGP2 = Double.parseDouble(lineArray[19]);
                    P_VSA_ppp_D = Double.parseDouble(lineArray[20]);
                    P_VSA_ppp_L = Double.parseDouble(lineArray[21]);
                    SdssC = Double.parseDouble(lineArray[22]);
                    SpMax_B_m = Double.parseDouble(lineArray[23]);
                    SpMax_B_s = Double.parseDouble(lineArray[24]);
                    SpMax4_Bh_m = Double.parseDouble(lineArray[25]);
                    SpMin2_Bh_m = Double.parseDouble(lineArray[26]);
                    SpPosA_B_v = Double.parseDouble(lineArray[27]);
                    Wap = Double.parseDouble(lineArray[28]);
                    X5v = Double.parseDouble(lineArray[29]);
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }
    }


    // additional methods for descriptors calculation
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
            int nH = 0;
            for (IAtom con : m.getConnectedAtomsList(at)) {
                if (con.getSymbol().equalsIgnoreCase("H"))
                    nH++;
            }
            if (nH > 0)
                return 2.152;

            return 1.779;
        }

        if (s.equalsIgnoreCase("H")) {
            IAtom connAt = m.getConnectedAtomsList(at).get(0);
            if (connAt.getSymbol().equalsIgnoreCase("O"))
                return 0.8;
            if (connAt.getSymbol().equalsIgnoreCase("N"))
                return 0.7;
            if (connAt.getSymbol().equalsIgnoreCase("P"))
                return 0.7;
            return 1.485;
        }

        return Descriptor.MISSING_VALUE;
    }

    private final static String[] PPP_TYPES = {
            "D",
            "A",
            "P",
            "N",
            "L",
            "Cyc"
    };

    private int Calculate_PRX() {

        int XAYbond = 0;

        for (int i=0; i<nAtoms; i++) {

            int nNbond = 0, nObond = 0, nOdbl = 0, nN_VD =0;

            String sym = CurMol.getAtom(i).getSymbol();

            if ( (sym.equalsIgnoreCase("C")) || (sym.equalsIgnoreCase("S")) ||
                    (sym.equalsIgnoreCase("P")) ) {

                for (int j=0; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j] != 0) {

                        if (ConnMatrix[j][j] == 7) {
                            // N
                            nNbond++;
                            for (int k=0; k<nAtoms; k++) {
                                if (k==j) continue;
                                if (ConnMatrix[j][k] != 0) nN_VD++;
                            }

                            // counts also H in N vd?
                            int nH;
                            try {
                                nH = CurMol.getAtom(j).getImplicitHydrogenCount();
                            } catch (Exception e) {
                                nH = 0;
                            }
                            nN_VD += nH;

                        }

                        if (ConnMatrix[j][j] == 8) {
                            // O
                            nObond++;
                            if (ConnMatrix[i][j] == 2) nOdbl++;
                        }
                    }
                }

                int nNObond = nObond + nNbond;

                // X-A-Y
                if ( nNObond == 2) {
                    XAYbond += 2;

                    // Corrections

                    // O-C-O
                    if ( (nObond == 2) && (nOdbl == 0) ) XAYbond--;

                    // >N-C-N<
                    if ( (nNbond == 2) && (nN_VD == 6) ) XAYbond--;

                    // >N-C-O or >N-C=O
                    if ( (nNbond == 1) && (nN_VD == 3) && (nObond == 1) ) XAYbond--;

                }

                if ( nNObond == 3) {
                    XAYbond += 3;
                    // Correction for Sulfonamide S(=O)(=O)N and Carboxamide
                    if ( (nOdbl > 0) && (((nNbond == 1) && (nN_VD == 3)) || ((nNbond == 2) && (nN_VD > 4))) ) XAYbond--;
                }

                if (nNObond == 4) {
                    XAYbond += 6;
                    if ( (sym.equalsIgnoreCase("S")) && (nOdbl > 0) &&
                            (((nNbond == 1) && (nN_VD == 3)) || ((nNbond == 2) && (nN_VD > 4)) || ((nNbond == 3)
                                    && (nN_VD > 6))) )
                        XAYbond--;
                }
            }
        }

        // NO
        int nNObond = 0;
        for (int i=0; i<nAtoms-1; i++) {
            if ( (ConnMatrix[i][i] == 7) || (ConnMatrix[i][i] == 8) )
                for (int j = i+1; j < nAtoms; j++)
                    if ( (ConnMatrix[j][j] == 7) || (ConnMatrix[j][j] == 8) ) {
                        if (ConnMatrix[i][j] != 0)
                            nNObond++;
                    }
        }

        return nNObond * 2 + XAYbond;
    }


    private boolean Check_HB() {

        for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
            IRing curRing = (IRing) MolRings.getAtomContainer(r);

            // consider aromatic PSA rings as in Dragon
            // Dragon seems NOT to consider if the ring structure is aromatic
            boolean AromaticPSA = true;
            for (int i = 0; i < curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));
                int VD = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (atIdx == j) continue;
                    if (ConnMatrix[atIdx][j] != 0)
                        VD++;
                }
                int nH;
                try {
                    nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }

                VD += nH;
                if (VD > 3) {
                    AromaticPSA = false;
                    break;
                }
            }

            if (!AromaticPSA)
                continue;

            // check all ring atoms for substituent
            // to find NH2 - OH in orto position
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));

                int VD = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (atIdx == j) continue;
                    if (ConnMatrix[atIdx][j] != 0)
                        VD++;
                }
                int nH;
                try {
                    nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }
                VD += nH;

                if (VD == 3) {

                    int FoundNH2 = Find_Attached_NH2(atIdx);
                    int FoundOH = Find_Attached_OH(atIdx);

                    // found a NH2 or OH group on atIdx
                    if ( (FoundOH > -1) || (FoundNH2 > -1) ) {

                        for (int atIdxB = 0; atIdxB < nAtoms; atIdxB++) {
                            if (atIdx == atIdxB) continue;

                            // check connected atoms in the same ring
                            if ( (ConnMatrix[atIdx][atIdxB] != 0) && (curRing.contains(CurMol.getAtom(atIdxB))) )  {

                                int FoundNH2_b = Find_Attached_NH2(atIdxB);
                                int FoundOH_b = Find_Attached_OH(atIdxB);

                                // check C=O not belonging to the same cycle
                                int FoundCO = -1;
                                for (int atIdxC = 0; atIdxC < nAtoms; atIdxC++) {
                                    if (atIdxC == atIdxB) continue;
                                    if (atIdxC == atIdx) continue;
                                    if (ConnMatrix[atIdxC][atIdxB] != 0) {

                                        if ( (!CurMol.getAtom(atIdxC).getFlag(CDKConstants.ISAROMATIC)) &&
                                                (!curRing.contains(CurMol.getAtom(atIdxC))) &&
                                                (ConnMatrix[atIdxC][atIdxC] == 6) ){
                                            int curVD = 0;
                                            int Odbl = 0;
                                            for (int k = 0; k < nAtoms; k++) {
                                                if (k == atIdxC) continue;
                                                if (ConnMatrix[atIdxC][k] != 0) {
                                                    curVD++;
                                                    if ( (ConnMatrix[atIdxC][k] == 2) && ConnMatrix[k][k] == 8)
                                                        Odbl++;
                                                }
                                            }
                                            if ( (curVD == 3) && (Odbl == 1) )
                                                FoundCO = atIdxC;
                                        }
                                    }
                                }

                                if ( (FoundOH > -1) && (FoundNH2_b > -1) )
                                    return true;
                                if ( (FoundOH > -1) && (FoundCO > -1) )
                                    return true;
                                if ( (FoundNH2 > -1) && (FoundCO > -1) )
                                    return true;
                            }
                        }
                    }
                }
            }

            // search NH2 - OH on different fused rings, and quinoline-type OH/NH2 -N=
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));

                int FoundNH2 = Find_Attached_NH2(atIdx);
                int FoundOH = Find_Attached_OH(atIdx);

                if ( (FoundNH2 > -1) || (FoundOH > -1)) {

                    // cycles on connected atoms to find one belonging to 2 different rings
                    for (int atIdxB = 0; atIdxB < nAtoms; atIdxB++) {
                        if (atIdx == atIdxB) continue;
                        if (ConnMatrix[atIdx][atIdxB] != 0) {

                            int countRings = 0;
                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                if (MolRings.getAtomContainer(rr).contains(CurMol.getAtom(atIdxB)))
                                    countRings++;
                            }
                            if (countRings != 2)
                                continue;

                            for (int atIdxC = 0; atIdxC < nAtoms; atIdxC++) {
                                if (atIdxC == atIdxB) continue;
                                if (atIdxC == atIdx) continue;
                                if ( (ConnMatrix[atIdxC][atIdxB] != 0) && (!curRing.contains(CurMol.getAtom(atIdxC)))) {
                                    // connected atoms at topo distance = 2 from atIdx, in a different but fused ring

                                    int FoundNH2_b = Find_Attached_NH2(atIdxC);
                                    int FoundOH_b = Find_Attached_OH(atIdxC);

                                    if ( (FoundNH2 > -1) && (FoundOH_b > -1) )
                                        return true;
                                    if ( (FoundOH > -1) && (FoundNH2_b > -1) )
                                        return true;

                                    // check quinolines
                                    if (ConnMatrix[atIdxC][atIdxC] == 7) {
                                        int nVD = 0, nSingle = 0, nDbl = 0;
                                        for (int k = 0; k < nAtoms; k++) {
                                            if (k == atIdxC) continue;
                                            if (ConnMatrix[atIdxC][k] != 0) {
                                                nVD++;
                                                if (ConnMatrix[atIdxC][k] == 1) nSingle++;
                                                if (ConnMatrix[atIdxC][k] == 2) nDbl++;
                                            }
                                        }
                                        if ( (nVD==2) && (nSingle==1) && (nDbl==2)) {
                                            // check ring sizes
                                            int RingSize1 = curRing.getAtomCount();
                                            int RingSize2 = 0;
                                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                                if (MolRings.getAtomContainer(rr).contains(CurMol.getAtom(atIdxC))) {
                                                    RingSize2 = MolRings.getAtomContainer(rr).getAtomCount();
                                                    break;
                                                }
                                            }
                                            if ( (RingSize1 == 6) && (RingSize2 == 6) )
                                                return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return false;
    }


    // Return atom index of attached NH2, or -1 if not found
    private int Find_Attached_NH2(int atIdx) {

        for (int i = 0; i < nAtoms; i++) {
            if (atIdx == i) continue;
            if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) continue;
            if ((ConnMatrix[atIdx][i] != 0)) {
                if (ConnMatrix[i][i] == 7) {
                    int nH;
                    try {
                        nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                    } catch (Exception e) {
                        nH = 0;
                    }
                    if (nH == 2)
                        return i;
                }
            }
        }

        return -1;
    }


    // Return atom index of attached OH, or -1 if not found
    private int Find_Attached_OH(int atIdx) {

        for (int i = 0; i < nAtoms; i++) {
            if (atIdx == i) continue;
            if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) continue;
            if ((ConnMatrix[atIdx][i] != 0)) {
                if (ConnMatrix[i][i] == 8) {
                    int nH;
                    try {
                        nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                    } catch (Exception e) {
                        nH = 0;
                    }
                    if (nH == 1)
                        return i;
                }
            }
        }

        return -1;
    }


    private double Check_NCS(IAtomContainer mol, int atom) {

        // Starting from S atom

        for (int i=0; i<nAtoms; i++) {
            if (i==atom) continue;
            if ( (ConnMatrix[i][atom]>0) && (AreEqual(GetName(mol,i),"C")) ) {
                for (int j=0; j<nAtoms; j++) {
                    if (j==i) continue;
                    if (AreEqual(GetName(mol,j),"N")) {
                        if (ConnMatrix[i][j]==2)
                            return 1;       // -N=C=S
                        if (ConnMatrix[i][j]==3)
                            return 0.5;     // N#C-S-
                    }
                }
            }
        }

        return 0;
    }


    private double Check_NCS() {

        double NCS = 0;

        for (int i=0; i<nAtoms; i++) {

            // Starting from C atom
            if (ConnMatrix[i][i] != 6) continue;

            int VD = 0;
            double bondS = 0;
            double bondN = 0;
            for (int j=0; j<nAtoms; j++) {
                if (i==j) continue;
                if (ConnMatrix[i][j] != 0) {
                    VD++;
                    if (ConnMatrix[i][j] == 7)
                        bondN = ConnMatrix[i][j];
                    if (ConnMatrix[i][j] == 16)
                        bondS = ConnMatrix[i][j];
                }
            }

            if (VD==2) {
                if ( (bondS == 2) && (bondN == 2) )
                    NCS += 1;
                if ( (bondS == 1) && (bondN == 3) )
                    NCS += 0.5;
            }

        }

        return NCS;
    }



    private int Check_BLM(IAtomContainer mol, int atom) {

        // Starting from N atom

        if (MolRings.contains(mol.getAtom(atom))) {
//            IRing R = (IRing)MolRings.getRings(mol.getAtom(atom)).getAtomContainer(0);
            for (IAtomContainer ac : MolRings.getRings(mol.getAtom(atom)).atomContainers() )  {
                IRing R = (IRing)ac;
                if (R.getAtomCount()==4) {
                    int N=0, C=0, CdblO=0, Other=0;
                    for (int i=0; i<4; i++) {
                        if (AreEqual(R.getAtom(i).getSymbol(),"N"))
                            N++;
                        else if (AreEqual(R.getAtom(i).getSymbol(),"C")) {
                            C++;
                            for (int k=0; k<nAtoms; k++)
                                if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[mol.indexOf(R.getAtom(i))][k]==2) ) {
                                    CdblO++;
                                }
                        } else
                            Other++;
                    }
                    if ( (N==1) && (C==3) && (CdblO==1) )
                        return 1;
                }
            }

        }

        return 0;
    }



    private int Calculate_POL() {

        int nPolSub = 0;

        for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
            IRing curRing = (IRing) MolRings.getAtomContainer(r);

            if (!IsRingPSA(curRing))
                continue;

            // check all ring atoms for substituent
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));
                if (ConnMatrix[atIdx][atIdx] == 6) {

                    for (int j=0; j<nAtoms; j++) {
                        if (atIdx == j) continue;
                        if ( (ConnMatrix[atIdx][j] != 0) ) {

                            // check if in same ring or fused ring
                            boolean sameRing = false;
                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                IRing bufRing = (IRing) MolRings.getAtomContainer(rr);
                                // Dragon seems not to consider fused rings if they are not PSA
                                if (!IsRingPSA(bufRing))
                                    continue;
                                if  ( (bufRing.contains(CurMol.getAtom(j))) && (bufRing.contains(CurMol.getAtom(atIdx))) ) {
                                    sameRing = true;
                                    break;
                                }
                            }
                            if (sameRing)
                                continue;

                            if (ConnMatrix[j][j] != 6) {
                                // all hetero atoms considered pol sub
                                nPolSub++;
                            } else {
                                int nHet = 0;
                                int cVD = 0;
                                for (int k=0; k<nAtoms; k++) {
                                    if (k==j) continue;
                                    if (ConnMatrix[k][j] != 0) {
                                        cVD++;
                                        if (ConnMatrix[k][k] != 6)
                                            nHet++;
                                    }
                                }

                                int nH;
                                try {
                                    nH = CurMol.getAtom(j).getImplicitHydrogenCount();
                                } catch (Exception e) {
                                    nH = 0;
                                }
                                cVD += nH;

                                if ( (cVD == 2) && (nHet >= 1) )
                                    nPolSub++;
                                if ( (cVD == 3) && (nHet >= 1) )
                                    nPolSub++;
                                if ( (cVD == 4) && (nHet >= 2) )
                                    nPolSub++;
                            }
                        }
                    }
                }
            }

        }

        return nPolSub;
    }


    private boolean IsRingPSA(IRing curRing) {

        // consider aromatic PSA rings as in Dragon
        // Dragon seems NOT to consider if the ring structure is aromatic
        boolean AromaticPSA = true;
        for (int i=0; i<curRing.getAtomCount(); i++) {
            int atIdx = CurMol.indexOf(curRing.getAtom(i));
//                if (!CurMol.getAtom(atIdx).getFlag(CDKConstants.ISAROMATIC)) {
//                    AromaticPSA = false; break;
//                }
            int VD = 0;
            for (int j=0; j<nAtoms; j++) {
                if (atIdx == j) continue;
                if (ConnMatrix[atIdx][j] != 0)
                    VD++;
            }
            int nH;
            try {
                nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
            } catch (Exception e) {
                nH = 0;
            }

            VD+= nH;
            if (VD > 3) {
                AromaticPSA = false;
                break;
            }
        }
        return AromaticPSA;
    }


    private int Calculate_RNG(IAtomContainer mol) {

        for (int i=0; i<MolRings.getAtomContainerCount(); i++) {

            IRing R = (IRing)MolRings.getAtomContainer(i);
            boolean SkipRing = false;

            // Skips benzene
            if (R.getAtomCount()==6)
                if (MoleculeUtilities.IsRingAromatic(R)) {
                    SkipRing = true;
                    for (int j=0; j<R.getAtomCount(); j++)
                        if (!(AreEqual(R.getAtom(j).getSymbol(),"C")))
                            SkipRing = false;
                }

            //
            if (!IsRingPSA(R))
                return 1;

            // Skips rings fused with benzene
            if (!SkipRing) {
                for (int at=0; at<R.getAtomCount(); at++) {
                    IAtom At = R.getAtom(at);
                    for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
                        if (r == i) continue;
                        IRing curR = (IRing)MolRings.getAtomContainer(r);
                        if (curR.contains(At)) {
                            // Checks if this ring is benzene
                            if ( (curR.getAtomCount() ==6 ) &&
                                    (MoleculeUtilities.IsRingAromatic(curR)) ) {
                                SkipRing = true;
                                break;
                            }
                        }
                    }
                    if (SkipRing) break;
                }
            }

            if (!SkipRing)
                return 1;
        }

        return 0;
    }



    private int Calculate_ALK(IAtomContainer mol) {

        boolean OnlyC;
        int nC=0;

        // Check if it is a hydrocarbon compound
        OnlyC = true;
        for (int i=0; i<nAtoms; i++) {
            if  (ConnMatrix[i][i] == 6) nC++;
            else OnlyC = false;
        }

        // Hydrocarbon compound
        if (OnlyC) {
            boolean HasAromaticBonds = false;
            for (int i=0; i<nAtoms; i++)
                for (int j=0; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j] == 1.5) {
                        HasAromaticBonds = true;
                        break;
                    }
                }

            if (HasAromaticBonds)
                return 0;

            // Checks alkane and alkene
            int nSingle=0, nDouble=0, nBonds=0;
            for (int i=0; i<(nAtoms-1); i++)
                for (int j=i; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j]>0) {
                        nBonds++;
                        if (ConnMatrix[i][j] == 1) nSingle++;
                        if (ConnMatrix[i][j] == 2) nDouble++;
                    }
                }
            if (nBonds == nSingle)
                return 1; // alkane

        }

        // Checks for hydrocarbon chains with at least 7 C
        if (nC>6) {

            for (int i=0; i<nAtoms; i++)
                if  (AreEqual(GetName(mol,i),"C")) {

                    IAtom chainAt = mol.getAtom(i);
                    int nH=0;
                    if (chainAt.getImplicitHydrogenCount()!=null)
                        nH = chainAt.getImplicitHydrogenCount();
                    if (nH!=3)
                        continue; // shall start from terminal CH3 atom

                    for (int j=0; j<nAtoms; j++)
                        if ((AreEqual(GetName(mol,j),"C"))&& (TopoMatrix[i][j]>6)) {

                            boolean ChainFound=true;

                            ShortestPaths shortestPaths = new ShortestPaths(mol, mol.getAtom(i));
                            List<IAtom> Path = Arrays.asList(shortestPaths.atomsTo(mol.getAtom(j)));

                            // DEPRECATED METHOD
//                            List<IAtom> Path = PathTools.getShortestPath(mol, mol.getAtom(i), mol.getAtom(j));

                            for (int k=1; k<Path.size(); k++) {
                                IAtom CurAtom = Path.get(k);
                                nH=0;
                                if (CurAtom.getImplicitHydrogenCount()!=null)
                                    nH = CurAtom.getImplicitHydrogenCount();
                                if ( (!(CurAtom.getSymbol().equalsIgnoreCase("C"))) ||
                                        (MolRings.contains(CurAtom)) ||
                                        (nH!=2) ) {
                                    // If atom is no Carbon or is not a chain, breaks
                                    ChainFound = false;
                                    break;
                                }
                            }
                            if (ChainFound)
                                return 1;
                        }
                }
        }


        return 0;
    }



    private int Calculate_AMP(IAtomContainer mol) {

        int count = 0;

        // Checks for COOH
        for (int i=0; i<nAtoms; i++) {

            boolean isC_COOH = false;
            if (!(AreEqual(GetName(mol,i),"C")))
                for (int j=0; j<nAtoms; j++) {
                    if ( (AreEqual(GetName(mol,j),"C")) && (ConnMatrix[i][j]==1) ) {
                        int Osingle=0, Odbl=0;
                        for (int k=0; k<nAtoms; k++) {
                            if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[j][k]==1) )
                                Osingle++;
                            if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[j][k]==2) )
                                Odbl++;
                        }
                        if ( (Osingle==1) && (Odbl==1) )
                            isC_COOH = true;
                    }
                }

            if (isC_COOH) {

                Boolean Pyridine=false, Amine=false;

                if (MolRings.contains(mol.getAtom(i))) {

                    // C atom is in ring, checks for Pyridine or aminobenzoic

                    boolean InvalidStruct = false;
                    IRing R = (IRing)MolRings.getRings(mol.getAtom(i)).getAtomContainer(0);
                    if ( (R.getAtomCount()==6) && (MoleculeUtilities.IsRingAromatic(R)) ) {
                        for (int j=0; j<6; j++) {
                            int Branches = 0;
                            int atomNum = mol.indexOf(R.getAtom(j));
                            if (atomNum!=i) {
                                if (AreEqual(GetName(mol,j),"C")) {
                                    for (int k=0; k<nAtoms; k++) {
                                        if ( (ConnMatrix[k][atomNum]>0) && (!(R.contains(mol.getAtom(k)))) ) {
                                            if ( (AreEqual(GetName(mol,k),"N")) && (GetHydrogens(mol, k)==2)) {
                                                Branches++;
                                                Amine = true;
                                            }
                                        }
                                    }
                                    if ((Branches>0)&&(!Amine)) {
                                        InvalidStruct = true;
                                        break;
                                    }
                                }
                                if (AreEqual(GetName(mol,j),"N")) {
                                    for (int k=0; k<nAtoms; k++) {
                                        if ( (ConnMatrix[k][atomNum]>0) && (!(R.contains(mol.getAtom(k)))) ) {
                                            InvalidStruct = true;
                                            break;
                                        }
                                    }
                                    Pyridine = true;
                                }
                            }
                        }

                        if ( (!InvalidStruct) && (Pyridine||Amine) )
                            count += 0.5;
                    }

                } else {

                    // No cyclic structure, checks for AA

                    int NH2=0, C=0, R=0;
                    for (int k=0; k<nAtoms; k++)
                        if (ConnMatrix[k][i]>0) {
                            if (AreEqual(GetName(mol,k),"N"))
                                if (GetHydrogens(mol, k)==2)
                                    NH2++;
                                else
                                    R++;
                            else if (AreEqual(GetName(mol,k),"C"))
                                C++;
                            else
                                R++;
                        }
                    if (NH2==1)
                        if ( ((C==1)&&(R==1)) || (C==2) )
                            count += 1;

                }

            }

        }

        return count;
    }



    ////// Utilities for this class

    private String GetName(IAtomContainer mol, int atom) {
        return mol.getAtom(atom).getSymbol();
    }

    private int GetHydrogens(IAtomContainer mol, int atom) {

        int nH=0;
        if (mol.getAtom(atom).getImplicitHydrogenCount()!=null)
            nH = mol.getAtom(atom).getImplicitHydrogenCount();

        return nH;
    }

    private boolean AreEqual(String AtomA, String AtomB) {
        return (AtomA.equalsIgnoreCase(AtomB));
    }

    private boolean isConjugated(IAtomContainer m, int atIdx, double[][] AdjConnAugMatrix) {

        int VD=0, nH=0, TotBondOrder=0;

        for (int i=0; i<m.getAtomCount(); i++) {
            if (i == atIdx) continue;
            if (AdjConnAugMatrix[i][atIdx] == 0) continue;

            VD++;
            TotBondOrder += m.getBond(m.getAtom(i), m.getAtom(atIdx)).getOrder().numeric();
        }
        try {
            nH = m.getAtom(atIdx).getImplicitHydrogenCount();
        } catch (Exception e) { }
        TotBondOrder += nH;
        VD += nH;

        boolean zero = false;
        int atEl = (int)AdjConnAugMatrix[atIdx][atIdx];
        if ( (atEl == 6) && (TotBondOrder == 4) && (VD == 3) )
            zero = true;
        if ( ((atEl == 7) || (atEl == 15)) && (TotBondOrder == 3) && (VD == 2) )
            zero = true;
        if ( ((atEl == 7) || (atEl == 15)) && (TotBondOrder == 4) && (VD == 3) )
            zero = true;
        if ( ((atEl == 5) || (atEl == 15)) && (TotBondOrder == 5) && (VD == 4) )
            zero = true;
        if ( (atEl == 5) && (TotBondOrder == 3) && (VD == 2) )
            zero = true;
        if ( (atEl == 16) && (TotBondOrder == 4) && (VD == 3) )
            zero = true;
        if ( ((atEl == 8) || (atEl == 16)) && (TotBondOrder == 2) && (VD == 1) )
            zero = true;

        if (zero) {

            for (int j=0; j<m.getAtomCount(); j++) {
                if (j == atIdx) continue;
                boolean first = false;
                if (AdjConnAugMatrix[j][atIdx] != 0) {

                    int curVD = 0, curnH = 0, curTotBondOrder = 0;
                    for (int k = 0; k < m.getAtomCount(); k++) {
                        if (k == j) continue;
                        if (AdjConnAugMatrix[j][k] != 0) {
                            curVD++;
                            curTotBondOrder += m.getBond(m.getAtom(k), m.getAtom(j)).getOrder().numeric();
                        }
                    }

                    try {
                        curnH = m.getAtom(j).getImplicitHydrogenCount();
                    } catch (Exception e) {
                    }
                    curTotBondOrder += curnH;
                    curVD += curnH;

                    int curAtEl = (int)AdjConnAugMatrix[j][j];

                    if ( (curAtEl == 6) && (curTotBondOrder == 4) && (curVD == 3) )
                        first = true;
                    if ( ((curAtEl == 7) || (curAtEl == 15)) && (curVD == 2) && (curTotBondOrder == 3) )
                        first = true;
                    if ( ((curAtEl == 7) || (curAtEl == 15)) && (curVD == 3) && (curTotBondOrder == 4) )
                        first = true;
                    if ( ((curAtEl == 5) || (curAtEl == 15)) && (curVD == 4) && (curTotBondOrder == 5) )
                        first = true;
                    if ( (curAtEl == 5) && (curVD == 2) && (curTotBondOrder == 3) )
                        first = true;
                    if ( (curAtEl == 16) && (curVD == 3) && (curTotBondOrder == 4) )
                        first = true;
                    if ( ((curAtEl == 8) || (curAtEl == 16)) && (curTotBondOrder == 2) && (curVD == 1) )
                        first = true;

                    if (first) {

                        for (int j2=0; j2<m.getAtomCount(); j2++) {
                            if (j2 == j) continue;
                            if (j2 == atIdx) continue;

                            if (AdjConnAugMatrix[j2][j] != 0) {

                                int cur2VD = 0, cur2nH = 0, cur2TotBondOrder = 0;
                                for (int kk = 0; kk < m.getAtomCount(); kk++) {
                                    if (kk == j2) continue;
                                    if (AdjConnAugMatrix[j2][kk] != 0) {
                                        cur2VD++;
//                                        cur2TotBondOrder += AdjConnAugMatrix[kk][j2];
                                        cur2TotBondOrder += m.getBond(m.getAtom(kk), m.getAtom(j2)).getOrder().numeric();
                                    }
                                }
                                try {
                                    cur2nH = m.getAtom(j2).getImplicitHydrogenCount();
                                } catch (Exception e) {
                                }
                                cur2TotBondOrder += cur2nH;
                                cur2VD += cur2nH;

                                int cur2AtEl = (int)AdjConnAugMatrix[j2][j2];

                                if ( (cur2AtEl == 6) && (cur2TotBondOrder == 4) && (cur2VD == 3) )
                                    return true;
                                if ( ((cur2AtEl == 7) || (cur2AtEl == 15)) && (cur2VD == 2) && (cur2TotBondOrder == 3) )
                                    return true;
                                if ( ((cur2AtEl == 7) || (cur2AtEl == 15)) && (cur2VD == 3) && (cur2TotBondOrder == 4) )
                                    return true;
                                if ( ((cur2AtEl == 5) || (cur2AtEl == 15)) && (cur2VD == 4) && (cur2TotBondOrder == 5) )
                                    return true;
                                if ( (cur2AtEl == 5) && (cur2VD == 2) && (cur2TotBondOrder == 3) )
                                    return true;
                                if ( (cur2AtEl == 16) && (cur2VD == 3) && (cur2TotBondOrder == 4) )
                                    return true;
                                if ( ((cur2AtEl == 8) || (cur2AtEl == 16)) && (cur2TotBondOrder == 2) && (cur2VD == 1) )
                                    return true;

                            }
                        }
                    }
                }
            }
        }

        return false;
    }

    /**
     * Sets CATS 2D atom types for each atom, as a list of string containing
     * all matching types for each atom.
     **/
    public ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        ArrayList[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false, tCyc=false;

            // Definition of CATS types
            //
            // A: O, N without H
            // N: [+], NH2
            // P: [-], COOH, POOH, SOOH

            // Hydrogens
            int H = 0;
            try {
                H = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("Unable to count H");
            }

            // counters
            int nSglBnd = 0, nOtherBnd = 0, VD = 0;
            int nC = 0, nDblO = 0, nOtherNonOBond=0, nSglOH = 0;
            for (int j=0; j<nSK; j++) {
                if (j==i) continue;
                if (ConnAugMatrix[i][j]>0) {

                    VD++;

                    if (ConnAugMatrix[j][j] == 6)
                        nC++;

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
                        nOtherBnd++;
                        if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
                            nDblO++;
                        else
                            nOtherNonOBond++;
                    }
                }

            }


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

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (H == 1) || (H ==2) )
                    tD = true;

            }

            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                    (CurAt.getFormalCharge() == 0) )  {

                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherNonOBond == 0) )
                    tN = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                if ( VD == nC)
                    tL = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                if ( (nC == 2) && ( (VD+H) == 2 ))
                    tL = true;
            }

            if (CurAt.isInRing())
                tCyc = true;


            // Sets final types
            if (tA) AtomTypes[i].add("A");
            if (tN) AtomTypes[i].add("N");
            if (tP) AtomTypes[i].add("P");
            if (tD) AtomTypes[i].add("D");
            if (tL) AtomTypes[i].add("L");
            if (tCyc) AtomTypes[i].add("Cyc");

        }

        return AtomTypes;
    }


}
