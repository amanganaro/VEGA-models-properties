package insilico.pparg_up.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.descriptor.blocks.P_VSA;
import insilico.descriptor.localization.StringSelectorDescriptors;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.smarts.SmartsPattern;
import insilico.core.descriptor.blocks.weights.iWeight;

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

    private final double MISSING_VALUE = -999;

    private IAtomContainer CurMol;
    private int nAtoms;
    private double[][] ConnMatrix;
    private int[][] TopoMatrix;
    private RingSet MolRings;

    public double ExperimentalValue;

    public double SpMAD_B_m = MISSING_VALUE;
    public double rGes = MISSING_VALUE;
    public double P_VSA_v_3 = MISSING_VALUE;
    public double P_VSA_LogP_4 = MISSING_VALUE;
    public double O = MISSING_VALUE;
    public double O_057 = MISSING_VALUE;
    public double nCconj = MISSING_VALUE;
    public double nCb = MISSING_VALUE;
    public double MLOGP2 = MISSING_VALUE;
    public double MLOGP = MISSING_VALUE;
    public double Mi = MISSING_VALUE;
    public double GATS8s = MISSING_VALUE;
    public double GATS1p = MISSING_VALUE;
    public double GATS1m = MISSING_VALUE;
    public double F06_C_O = MISSING_VALUE;
    public double D_Dtr06 = MISSING_VALUE;
    public double CATS2D_07_AL = MISSING_VALUE;
    public double CATS2D_03_LL = MISSING_VALUE;
    public double C = MISSING_VALUE;
    public double C_026 = MISSING_VALUE;
    public double B07_C_O = MISSING_VALUE;
    public double ATSC7m = MISSING_VALUE;

    public double[] getDescriptors() {
        return new double[]{SpMAD_B_m, rGes, P_VSA_v_3, P_VSA_LogP_4, O, O_057, nCconj, nCb, MLOGP2, MLOGP,
                Mi, GATS8s, GATS1p, GATS1m, F06_C_O, D_Dtr06, CATS2D_07_AL, CATS2D_03_LL, C, C_026, B07_C_O, ATSC7m
        };
    }

    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        if (fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateSPMad(mol);
        CalculatePVSA(mol);
        CalculateOPerc(mol);
        CalculateACF(mol);
        CalculateFunctionalGroups(mol);
        CalculateMLogP(mol);
        CalculateMi(mol);
        CalculateAutoCorrelation(mol);
        CalculateGats8S(mol);

    }

    private void CalculateGats8S(InsilicoMolecule mol) {
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();

        WeightsIState curWeight = new WeightsIState();
        double[] w = curWeight.getWeights(m, true);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        int lag = 8;

        double GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag) {

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


        GATS8s = GearyAC;
        System.out.println();
    }

    private void CalculateAutoCorrelation(InsilicoMolecule mol) {

//        GATS1m; GATS1p; GATS8s;

        try {
            IAtomContainer m;
            try {
                IAtomContainer orig_m = mol.GetStructure();
                m = Manipulator.AddHydrogens(orig_m);
            } catch (InvalidMoleculeException | GenericFailureException e) {
                log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
                return;
            }

            // Gets matrices
            int[][] TopoMatrix;
            try {
                TopoMatrix = TopoDistanceMatrix.getMatrix(m);
            } catch (Exception e) {
                log.warn(e.getMessage());
                return;
            }

            int nSK = m.getAtomCount();

            // Cycle for all found weighting schemes
            iBasicWeight curWeight = new WeightsMass();
            double[] w = curWeight.getScaledWeights(m);

            // Calculates weights averages
            double wA = 0;
            for (int i=0; i<nSK; i++)
                wA += w[i];
            wA = wA / ((double) nSK);

//            for (int lag=1; lag<=8; lag++) {

            int lag = 1;
            double GearyAC = 0;
            double denom = 0, delta = 0;

            for (int i = 0; i < nSK; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j = 0; j < nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        GearyAC += Math.pow((w[i] - w[j]), 2);
                        delta++;
                    }
            }

            if (delta > 0) {
                if (denom == 0) {
                    GearyAC = 0;
                } else {
                    GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double) (nSK - 1))) * denom);
                }
            }

            GATS1m = GearyAC;


            curWeight = new WeightsPolarizability();
            w = curWeight.getScaledWeights(m);

            // Calculates weights averages
            wA = 0;
            for (int i=0; i<nSK; i++)
                wA += w[i];
            wA = wA / ((double) nSK);

//            for (int lag=1; lag<=8; lag++) {

            lag = 1;
            GearyAC = 0;
            denom = 0;
            delta = 0;

            for (int i = 0; i < nSK; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j = 0; j < nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        GearyAC += Math.pow((w[i] - w[j]), 2);
                        delta++;
                    }
            }

            if (delta > 0) {
                if (denom == 0) {
                    GearyAC = 0;
                } else {
                    GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double) (nSK - 1))) * denom);
                }
            }


            GATS1p = GearyAC;
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateMi(InsilicoMolecule mol) {

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(String.format(StringSelectorCore.getString("descriptors_invalid_structure"), mol.GetSMILES()));
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nTotH=0;



            //// Counts on atoms

            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];
            }


            for (int w=0; w<5; w++) {

                iBasicWeight ws;
                switch (w) {
                    case 0:
                        ws = new WeightsMass();
                        break;
                    case 1:
                        ws = new WeightsVanDerWaals();
                        break;
                    case 2:
                        ws = new WeightsPolarizability();
                        break;
                    case 3:
                        ws = new WeightsElectronegativity();
                        break;
                    case 4:
                        ws = new WeightsIonizationPotential();
                        break;
                    default:
                        throw new Exception(StringSelectorCore.getString("descriptors_weight_not_found"));
                }

                double[] weights = ws.getScaledWeights(curMol);
                double weightH = ws.getScaledWeight("H");

                double sum = 0;
                for (int i=0; i<nSK; i++) {
                    if (weights[i] == Descriptor.MISSING_VALUE) {
                        sum = Descriptor.MISSING_VALUE;
                        break;
                    } else {
                        // all values INCLUDING MW are scaled on carbon
                        sum += weights[i];
                        if (H[i] > 0)
                            sum += weightH * H[i];
                    }
                }

                double ave = Descriptor.MISSING_VALUE;
                if (sum != Descriptor.MISSING_VALUE)
                    ave = sum/(nSK + nTotH);

                switch (w) {
                    case 0:
                    case 1:
                    case 2:
                    case 3:
                        break;
                    case 4:
                        Mi = ave;
                        break;
                    default:
                        throw new Exception(StringSelectorCore.getString("descriptors_weight_not_found"));
                }
            }

        } catch (Throwable e) {
            log.warn(String.format(StringSelectorCore.getString("descriptors_unable_calculate"), mol.GetSMILES(), e.getMessage()));
        }
    }

    private void CalculateMLogP(InsilicoMolecule mol) {

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

        MLOGP = LogP;
        MLOGP2 = Math.pow(MLOGP, 2);
    }

    private void CalculateFunctionalGroups(InsilicoMolecule mol) {

        String[][] FG_SMARTS = {
                {"10", "nCb–", StringSelectorCore.getString("descriptors_fg_nCb"), "[$([c;D3]1ccccc1)]"},
        };

        try {
            int FGNumber = FG_SMARTS.length;
            SmartsPattern[] Queries = new SmartsPattern[FGNumber];
            extractDescriptor(FG_SMARTS, FGNumber, Queries);

            boolean status = Queries[0].matches(mol.GetStructure());
            List<Mappings> mappings = new ArrayList<>();
            int nmatch = 0;
            if (status) {
                mappings.add(Queries[0].matchAll(mol.GetStructure()));
                nmatch = Queries[0].matchAll(mol.GetStructure()).countUnique();
            }

            nCb = nmatch;

            FG_SMARTS = new String[][]{
                    {"11", "nCconj", StringSelectorCore.getString("descriptors_fg_nCconj"), "[$(C=CC=*),$(C(=*)C=*),$(C(=*)[a]),$(C=[*][a])]"},
            };

            FGNumber = FG_SMARTS.length;
            Queries = new SmartsPattern[FGNumber];
            extractDescriptor(FG_SMARTS, FGNumber, Queries);


            status = Queries[0].matches(mol.GetStructure());
            mappings = new ArrayList<>();
            nmatch = 0;
            if (status) {
                mappings.add(Queries[0].matchAll(mol.GetStructure()));
                nmatch = Queries[0].matchAll(mol.GetStructure()).countUnique();
            }

            nCconj = nmatch;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


    }

    private void CalculateACF(InsilicoMolecule mol) {
        IAtomContainer CurMol;
        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        // Gets fragments
        GhoseCrippenACF GC = new GhoseCrippenACF(CurMol, false);
        int[] Frags = GC.GetACF();

        // Count fragments to fill descriptors
        int[] FragDescriptors = new int[GhoseCrippenACF.ACF_NAMES.length];
        for (int i : FragDescriptors)
            i = 0;
        for (int frag : Frags) {
            if (frag == Descriptor.MISSING_VALUE)
                continue;
            FragDescriptors[frag]++;
        }

        // Count fragments not mapped directly into atoms (H)
        for (int key : GC.getNotMappedFragCount().keySet())
            FragDescriptors[key] = GC.getNotMappedFragCount().get(key);

        // Set descriptors
        for (int i = 0; i < GhoseCrippenACF.ACF_NAMES.length; i++) {
            if (i == 26)
                C_026 = FragDescriptors[i];
            if(i == 57)
                O_057 = FragDescriptors[i];
        }
    }

    private void CalculateOPerc(InsilicoMolecule mol) {


        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(String.format(StringSelectorCore.getString("descriptors_invalid_structure"), mol.GetSMILES()));
            return;
        }

        try {
            int nSK = curMol.getAtomCount();
            int[] H = new int[nSK];

            int nTotH=0;
            int nO=0;


            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];


                if (CurAt.getSymbol().equalsIgnoreCase("O"))
                    nO++;

            }

            O = (nO/(double)(nSK + nTotH))*100;
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


    }

    private void CalculatePVSA(InsilicoMolecule mol) {
        DescriptorBlock block = new P_VSA();
        try {
            block.Calculate(mol);
            P_VSA_LogP_4 = block.GetByName("P_VSA_LogP_4").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        //p vsa v 3
        IAtomContainer m = null;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
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

        double[] w = (new WeightsVanDerWaals()).getScaledWeights(m);
        for (int b=1; b<=4; b++) {
            double PVSA = 0;
            for (int i = 0; i < nSK; i++) {
                if (w[i] == Descriptor.MISSING_VALUE) continue;
                if (b == CalculateBin("v", w[i]))
                    PVSA += VSA[i];
            }
            if(b==3)
                P_VSA_v_3 = PVSA;
        }
    }

    private void CalculateSPMad(InsilicoMolecule mol) {
        try {
            IAtomContainer curMol;
            try {
                curMol = mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
                return;
            }

            // Adj matrix available for calculations
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = mol.GetMatrixBurden();
            double[] w = (new WeightsMass()).getScaledWeights(curMol);

            for (int i=0; i<nSK; i++)
                Mat[i][i] = w[i];

            Matrix DataMatrix = new Matrix(Mat);
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            double[] eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);

            double EigAve = 0;
            double EigMax = eigenvalues[0];
            double EigMin = eigenvalues[0];
            for (double val : eigenvalues) {
                EigAve += val;
                if (val > EigMax)
                    EigMax = val;
                if (val< EigMin)
                    EigMin = val;
            }
            EigAve = EigAve / (double) nSK;

            double EigDev = 0;
            for (double val : eigenvalues)
                EigDev += Math.abs(val - EigAve);
            SpMAD_B_m = EigDev / (double) nSK;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

    }

    public void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {
        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-PPARg_up/src/main/resources/data/dataset.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            br.readLine();
            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[2].trim()).GetSMILES())) {

                    ExperimentalValue = Double.parseDouble(lineArray[14]);

                    SpMAD_B_m = Double.parseDouble(lineArray[15]);
                    rGes = Double.parseDouble(lineArray[16]);
                    P_VSA_v_3 = Double.parseDouble(lineArray[17]);
                    P_VSA_LogP_4 = Double.parseDouble(lineArray[18]);
                    O = Double.parseDouble(lineArray[19]);
                    O_057 = Double.parseDouble(lineArray[20]);
                    nCconj = Double.parseDouble(lineArray[21]);
                    nCb = Double.parseDouble(lineArray[22]);
                    MLOGP2 = Double.parseDouble(lineArray[23]);
                    MLOGP = Double.parseDouble(lineArray[24]);
                    Mi = Double.parseDouble(lineArray[25]);
                    GATS8s = Double.parseDouble(lineArray[26]);
                    GATS1p = Double.parseDouble(lineArray[27]);
                    GATS1m = Double.parseDouble(lineArray[28]);
                    F06_C_O = Double.parseDouble(lineArray[29]);
                    D_Dtr06 = Double.parseDouble(lineArray[30]);
                    CATS2D_07_AL = Double.parseDouble(lineArray[31]);
                    CATS2D_03_LL = Double.parseDouble(lineArray[32]);
                    C = Double.parseDouble(lineArray[33]);
                    C_026 = Double.parseDouble(lineArray[34]);
                    B07_C_O = Double.parseDouble(lineArray[35]);
                    ATSC7m = Double.parseDouble(lineArray[36]);
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }

    }

    private int CalculateBin(String curWeight, double value) {

        if (curWeight.equalsIgnoreCase("v")) {
            // VdW volume
            if (value < 0.5) return 1;
            if (value < 1.0) return 2;
            if (value < 1.3) return 3;
            return 4;
        }


        return 0;
    }


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

    private void extractDescriptor(String[][] FG_SMARTS, int FGNumber, SmartsPattern[] queries) {
        for (int i = 0; i < FGNumber; i++) {
            try {
                queries[i] = SmartsPattern.create(FG_SMARTS[i][3]).setPrepare(false);
            } catch (Exception e) {
                log.warn("Unable to parse SMARTS in functional groups " + FG_SMARTS[i][3]);
                queries[i] = null;
            }
        }
    }

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




}
