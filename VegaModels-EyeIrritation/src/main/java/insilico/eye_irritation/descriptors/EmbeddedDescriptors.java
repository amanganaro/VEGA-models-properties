package insilico.eye_irritation.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.Rings;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.iWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsHydrophobicityGC;
import insilico.core.descriptor.blocks.weights.other.WeightsMolarRefractivityGC;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.eye_irritation.descriptors.weights.WeightsIState;
//import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.*;

public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private final int MISSING_VALUE = -999;
    public double SpMaxA_AEA_dm = MISSING_VALUE;
    public double B04_N_O = MISSING_VALUE;
    public double MATS8s = MISSING_VALUE;
    public double P_VSA_e_5 = MISSING_VALUE;
    public double CATS2D_09_DA = MISSING_VALUE;
    public double N_067 = MISSING_VALUE;
    public double Eig07_AEA_dm = MISSING_VALUE;
    public double nOxiranes = MISSING_VALUE;
    public double CATS2D_02_DN = MISSING_VALUE;
    public double MATS3s = MISSING_VALUE;
    public double nR_Ct = MISSING_VALUE;
    public double nR_CH_X = MISSING_VALUE;
    public double CATS2D_06_PP = MISSING_VALUE;
    public double B05_N_S = MISSING_VALUE;
    public double N_068 = MISSING_VALUE;
    public double F06_N_F = MISSING_VALUE;
    public double B10_C_S = MISSING_VALUE;
    public double B08_N_S = MISSING_VALUE;
    public double ATS5i = MISSING_VALUE;
    public double ATSC1m = MISSING_VALUE;
    public double B06_C_N = MISSING_VALUE;
    public double nCH2RX = MISSING_VALUE;
    public double CATS2D_03_AP = MISSING_VALUE;
    public double CATS2D_04_DN = MISSING_VALUE;
    public double B04_C_O = MISSING_VALUE;
    public double nSO3OH = MISSING_VALUE;
    public double B01_N_S = MISSING_VALUE;
    public double GATS4p = MISSING_VALUE;
    public double MATS7m = MISSING_VALUE;
    public double B02_C_C = MISSING_VALUE;
    public double SpDiam_EA_dm = MISSING_VALUE;
    public double nRCOOR = MISSING_VALUE;
    public double B02_O_S = MISSING_VALUE;
    public double SpMaxA_EA_dm = MISSING_VALUE;
    public double B04_C_S = MISSING_VALUE;
    public double B10_C_N = MISSING_VALUE;
    public double P_VSA_m_2 = MISSING_VALUE;
    public double Eig05_EA_dm = MISSING_VALUE;
    public double F04_N_Br = MISSING_VALUE;
    public double B07_C_S = MISSING_VALUE;
    public double CATS2D_08_AN = MISSING_VALUE;
    public double Cl_088 = MISSING_VALUE;
    public double CATS2D_01_DN = MISSING_VALUE;
    public double B03_O_S = MISSING_VALUE;
    public double GATS5m = MISSING_VALUE;
    public double SM15_AEA_dm = MISSING_VALUE;
    public double SM03_AEA_dm = MISSING_VALUE;
    public double F05_N_N = MISSING_VALUE;
    public double B02_O_O = MISSING_VALUE;
    public double B01_C_Br = MISSING_VALUE;


    public EmbeddedDescriptors(InsilicoMolecule Mol) throws Exception {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws Exception {

        CalculateFunctionalGroup(Mol);
        CalculateAtomCenterFragmentsDescriptors(Mol);
        CalculatePVSADescriptors(Mol);
        CalculateCATS2DDescriptors(Mol);
        CalculateAutocorrelationDescriptors(Mol);
        CalculateEdgeAdiacencyDescriptors(Mol);
        CalculateAtomPairs2DDescriptors(Mol);

    }

    private void CalculateFunctionalGroup(InsilicoMolecule mol) throws Exception {
        String[][] FG_SMARTS = {
                {"nR=Ct", "[$([C;D3]([#6])([#6])=C)]"},
                {"nRCOOR", "O=[$([C;D2]),$([C;D3]C)][O;D2][C,c]"},
                {"nSO3OH", "[S,O;D1;!-]S([S,O])(=[S,O])=[S,O]"},
                {"nCH2RX", "[C,c][C;D2;!R][Cl,Br,F,I]"},
                {"nR=CHX", "C=[C;D2;!R][Cl,Br,F,I]"},
                {"nOxiranes", "C1CO1"},
        };
        // Init SMARTS
        int FGNumber = FG_SMARTS.length;
        Pattern[] Queries = new Pattern[FGNumber];
        for (int i = 0; i < FGNumber; i++) {
            try {

                String curSmarts = FG_SMARTS[i][1];

                Queries[i] = SmartsPattern.create(curSmarts).setPrepare(false);
            } catch (Exception e) {
                log.warn("descriptors_parsing_smarts_error");
                Queries[i] = null;
            }
        }

        for (int i=0; i<FG_SMARTS.length; i++) {

            int nmatch = 0;
            List<Mappings> mappings = new ArrayList<>();
            boolean status;

            try {
                status = Queries[i].matches(mol.GetStructure());
                if (status) {
                    mappings.add(Queries[i].matchAll(mol.GetStructure()));
                    nmatch = Queries[i].matchAll(mol.GetStructure()).countUnique();
                }
            } catch (Exception e) {
                log.warn(String.format(StringSelectorCore.getString("descriptors_fgquery_init_fail"), i+1, e.getMessage()));
                continue;
            }

            if (i==0)
                nR_Ct = nmatch;
            if (i==1)
                nRCOOR = nmatch;
            if (i==2)
                nSO3OH = nmatch;
            if (i==3)
                nCH2RX = nmatch;
            if (i==4)
                nR_CH_X = nmatch;
            if (i==5)
                nOxiranes = nmatch;

        }

    }
    private void CalculateAtomCenterFragmentsDescriptors(InsilicoMolecule mol) throws InvalidMoleculeException {
        IAtomContainer CurMol;
        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        // Gets fragments
        GCAtomCentredFragments GC = new GCAtomCentredFragments(CurMol, false);
        int[] Frags = GC.GetACF();

        // Count fragments to fill descriptors
        int[] FragDescriptors = new int[121];
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
        N_067 = FragDescriptors[67];
        N_068 = FragDescriptors[68];
        Cl_088 = FragDescriptors[88];

    }
    private void CalculatePVSADescriptors(InsilicoMolecule mol){

        String W_MASS = (new WeightsMass()).getSymbol();
        String W_ELECTRONEGATIVITY = (new WeightsElectronegativity()).getSymbol();

        String[] WEIGHTS = {
                W_MASS,
                W_ELECTRONEGATIVITY,
        };

        // P_VSA are calculated on H filled molecules
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("invalid_structure");
            return;
        }

        int nSK = m.getAtomCount();
        int nBO = m.getBondCount();

        GhoseCrippenACF GC = new GhoseCrippenACF(m, true);

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

            // TODO: cosa si usa al posto del metodo deprecato
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


        // Cycle for all found weighting schemes
        for (String curWeight : WEIGHTS) {

            // Sets needed weights
            double[] w = null;

            if (curWeight.equalsIgnoreCase(W_MASS)) {
                w = (new WeightsMass()).getScaledWeights(m);
                P_VSA_m_2 = 0;

                for (int i = 0; i < nSK; i++) {
                    if (w[i] == Descriptor.MISSING_VALUE) continue;
                    if (2 == CalculateBin(curWeight, w[i]))
                        P_VSA_m_2 += VSA[i];
                }
            } else if (curWeight.equalsIgnoreCase(W_ELECTRONEGATIVITY)) {
                w = (new WeightsElectronegativity()).getScaledWeights(m);
                P_VSA_e_5 = 0;
                for (int i = 0; i < nSK; i++) {
                    if (w[i] == Descriptor.MISSING_VALUE) continue;
                    if (5 == CalculateBin(curWeight, w[i]))
                        P_VSA_e_5 += VSA[i];
                }
            }
        }
    }
    private void CalculateCATS2DDescriptors(InsilicoMolecule mol){

        String TYPE_D = "D";
        String TYPE_A = "A";
        String TYPE_P = "P";
        String TYPE_N = "N";
        String TYPE_L = "L";

        String[][] AtomCouples = {
                {TYPE_D, TYPE_D},
                {TYPE_D, TYPE_A},
                {TYPE_D, TYPE_P},
                {TYPE_D, TYPE_N},
                {TYPE_D, TYPE_L},
                {TYPE_A, TYPE_A},
                {TYPE_A, TYPE_P},
                {TYPE_A, TYPE_N},
                {TYPE_A, TYPE_L},
                {TYPE_P, TYPE_P},
                {TYPE_P, TYPE_N},
                {TYPE_P, TYPE_L},
                {TYPE_N, TYPE_N},
                {TYPE_N, TYPE_L},
                {TYPE_L, TYPE_L},
        };

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("invalid_structure");
            return;
        }

        int nSK = curMol.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        double[][] ConnAugMatrix = null;
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(curMol, ConnAugMatrix);


        for (String[] atomCouple : AtomCouples) {

            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                for (int j = i; j < nSK; j++) {
                    if (isInCatsType(CatsTypes[i], CatsTypes[j], atomCouple[0], atomCouple[1]))
                        if (TopoMat[i][j] < 10)
                            desc[TopoMat[i][j]]++;
                }
            }

            for (int i = 1; i < desc.length; i++) {
                if ((i==1) && (atomCouple[0].equals("D")) && (atomCouple[1].equals("N")))
                    CATS2D_01_DN = desc[i];
                if ((i==2) && (atomCouple[0].equals("D")) && (atomCouple[1].equals("N")))
                    CATS2D_02_DN = desc[i];
                if ((i==3) && (atomCouple[0].equals("A")) && (atomCouple[1].equals("P")))
                    CATS2D_03_AP = desc[i];
                if ((i==4) && (atomCouple[0].equals("D")) && (atomCouple[1].equals("N")))
                    CATS2D_04_DN = desc[i];
                if ((i==6) && (atomCouple[0].equals("P")) && (atomCouple[1].equals("P")))
                    CATS2D_06_PP = desc[i];
                if ((i==8) && (atomCouple[0].equals("A")) && (atomCouple[1].equals("N")))
                    CATS2D_08_AN = desc[i];
                if ((i==9) && (atomCouple[0].equals("D")) && (atomCouple[1].equals("A")))
                    CATS2D_09_DA = desc[i];
            }
        }
    }

    private void CalculateAutocorrelationDescriptors(InsilicoMolecule mol){

        ArrayList<iWeight> bWeights = new ArrayList<>();
        bWeights.add(new WeightsMass());
        bWeights.add(new WeightsPolarizability());
        bWeights.add(new WeightsIonizationPotential());
        bWeights.add(new WeightsIState());

        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("invalid_structure");
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

          for (iWeight curWeight : bWeights) {
            double[] w;

            if (curWeight.getClass() == WeightsIState.class) {

                // I-States
                w = ((WeightsIState)curWeight).getWeights(m, true);

                // correction for compatibility with D7
                // H I-state is always 1
                for (int i=0; i<nSK; i++) {
                    if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                        w[i] = 1;
                }

            } else {

                // All other weights are basic weights (scaled values)
                w = ((iBasicWeight) curWeight).getScaledWeights(m);
            }

            // If one or more weights are not available, sets all to missing value
            boolean MissingWeight = false;
            for (int i=0; i<nSK; i++)
                if (w[i] == Descriptor.MISSING_VALUE) {
                    MissingWeight = true;
                    break;
                }
            if (MissingWeight)
                return;

            // Calculates weights averages
            double wA = 0;
            for (int i=0; i<nSK; i++)
                wA += w[i];
            wA = wA / ((double) nSK);

            // Calculates autocorrelations
            for (int lag=1; lag<=8; lag++) {

                double AC=0, ACS=0, MoranAC=0, GearyAC=0;
                double denom = 0, delta = 0;

                for (int i=0; i<nSK; i++) {

                    denom += Math.pow((w[i] - wA), 2);

                    for (int j=0; j<nSK; j++)
                        if (TopoMatrix[i][j] == lag) {
                            AC += w[i] * w[j];
                            ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                            MoranAC += (w[i] - wA) * (w[j] - wA);
                            GearyAC += Math.pow((w[i] - w[j]), 2);
                            delta++;
                        }
                }

                if (delta > 0) {
                    if (denom == 0) {
                        MoranAC = 1;
                        GearyAC = 0;
                    } else {
                        MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                        GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
                    }
                }

                // AC transformed in log form
                AC /= 2.0;
                AC = Math.log(1 + AC);

                ACS /= 2.0;

                // Sets descriptors
                if ((curWeight.getSymbol().equals("m")) && (lag == 1))
                    ATSC1m = ACS;

                if ((curWeight.getSymbol().equals("m")) && (lag == 5))
                    GATS5m = GearyAC;

                if ((curWeight.getSymbol().equals("m")) && (lag == 7))
                    MATS7m = MoranAC;

                if ((curWeight.getSymbol().equals("s")) && (lag == 3))
                    MATS3s = MoranAC;

                if ((curWeight.getSymbol().equals("s")) && (lag == 8))
                    MATS8s = MoranAC;

                if ((curWeight.getSymbol().equals("i")) && (lag == 5))
                    ATS5i = AC;

                if ((curWeight.getSymbol().equals("p")) && (lag == 4))
                    GATS4p = GearyAC;

            }
        }

    }
    private void CalculateAtomPairs2DDescriptors(InsilicoMolecule mol){
        int MAX_TOPO_DISTANCE = 10;
        String[] ATOMS_FOR_COUPLES = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si"};
        String[][] ATOM_COUPLES;

        ArrayList<String[]> bufList = new ArrayList<>();
        for (int i=0; i<ATOMS_FOR_COUPLES.length; i++)
            for (int j=i; j<ATOMS_FOR_COUPLES.length; j++) {
                String[] buf = new String[2];
                buf[0] = ATOMS_FOR_COUPLES[i];
                buf[1] = ATOMS_FOR_COUPLES[j];
                bufList.add(buf);
            }

        ATOM_COUPLES = new String[bufList.size()][2];
        for (int i=0; i<bufList.size(); i++) {
            ATOM_COUPLES[i][0] = bufList.get(i)[0];
            ATOM_COUPLES[i][1] = bufList.get(i)[1];
        }

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn("Invalid structure");
            return;
        }

        for (String[] atomCouple : ATOM_COUPLES) {

            int[] descB = new int[MAX_TOPO_DISTANCE];
            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCouple[1])) {

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j] - 1] = 1;
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


            for (int i = 0; i < descB.length; i++) {
                int lag = i + 1;

                if (lag == 1) {
                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B01_N_S = descB[i];
                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("Br"))))
                        B01_C_Br = descB[i];
                }

                if (lag == 2) {
                    if ((atomCouple[0].equalsIgnoreCase("O") && (atomCouple[1].equalsIgnoreCase("O"))))
                        B02_O_O = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("C"))))
                        B02_C_C = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("O") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B02_O_S = descB[i];
                }

                if (lag == 3) {
                    if ((atomCouple[0].equalsIgnoreCase("O") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B03_O_S = descB[i];
                }

                if (lag == 4) {
                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("O"))))
                        B04_C_O = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B04_C_S = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("O"))))
                        B04_N_O = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("Br"))))
                        F04_N_Br = descF[i];
                }

                if (lag == 5) {
                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B05_N_S = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("N"))))
                        F05_N_N = descF[i];
                }

                if (lag == 6) {
                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("N"))))
                        B06_C_N = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("F"))))
                        F06_N_F = descF[i];
                }

                if (lag == 7) {
                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B07_C_S = descB[i];
                }

                if (lag == 8) {
                    if ((atomCouple[0].equalsIgnoreCase("N") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B08_N_S = descB[i];
                }

                if (lag == 10) {
                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("S"))))
                        B10_C_S = descB[i];

                    if ((atomCouple[0].equalsIgnoreCase("C") && (atomCouple[1].equalsIgnoreCase("N"))))
                        B10_C_N = descF[i];
                }

            }


        }

    }
    private void CalculateEdgeAdiacencyDescriptors(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("invalid_structure");
            return;
        }

        int nSK = curMol.getAtomCount();
        int nBO = curMol.getBondCount();

        // Only for mol with nSK>1
        if (nSK < 2) {
            log.warn("less_2_atoms_error");
            return;
        }

        String[] MATRICES = new String[]{"EA", "AEA"};
        String[] NAMES = new String[]{"spMaxA","SM15", "SM03", "Eig07", "Eig05","SpDiam"};

        for (String curMat : MATRICES) {

            // Gets basic matrix
            double[][][] EdgeAdjMat = null;
            try {
                EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
            } catch (GenericFailureException e) {
                log.warn(e.getMessage());
                return;
            }

            double[][] curDataMatrix = new double[nBO][nBO];

            // plain EA
            for (int i = 0; i < nBO; i++)
                for (int j = 0; j < nBO; j++)
                    curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

            // weights for the matrix
            double[] w = new double[nBO];

            //// Dipole moment
            for (int i = 0; i < nBO; i++) {
                IAtom a = curMol.getBond(i).getAtom(0);
                IAtom b = curMol.getBond(i).getAtom(1);
                double CurVal = GetDipoleMoment(curMol, a, b);
                if (CurVal == 0)
                    CurVal = GetDipoleMoment(curMol, b, a);
                w[i] = CurVal;
            }

            if (curMat.equalsIgnoreCase("EA")) {
                for (int i = 0; i < nBO; i++) {
                    for (int j = 0; j < nBO; j++) {
                        if (curDataMatrix[i][j] != 0)
                            curDataMatrix[i][j] = w[j];
                    }
                }
            }

            // Build AEA-based matrix
            // Replace just diagonal EA elements with value of the bond weight
            if (curMat.equalsIgnoreCase("AEA")) {
                for (int i = 0; i < nBO; i++)
                    curDataMatrix[i][i] = w[i];
            }

            // Calculate standard eigenvalue-based descriptors
            EigenvalueBasedDescriptors eigDesc = new EigenvalueBasedDescriptors();
            double[] descriptors = eigDesc.Calculate(curDataMatrix, nBO);

            for (String curName : NAMES) {

                if (curName.equalsIgnoreCase("spMaxA") && curMat.equalsIgnoreCase("AEA"))
                    SpMaxA_AEA_dm = descriptors[0];

                if (curName.equalsIgnoreCase("spMaxA") && curMat.equalsIgnoreCase("EA"))
                    SpMaxA_EA_dm = descriptors[0];

                if (curName.equalsIgnoreCase("SM15") && (curMat.equalsIgnoreCase("AEA")))
                    SM15_AEA_dm = descriptors[1];

                if (curName.equalsIgnoreCase("SM03") && (curMat.equalsIgnoreCase("AEA")))
                    SM03_AEA_dm = descriptors[2];

                if (curName.equalsIgnoreCase("Eig07") && (curMat.equals("AEA")))
                    Eig07_AEA_dm = descriptors[3];

                if (curName.equalsIgnoreCase("Eig05") && (curMat.equals("EA")))
                    Eig05_EA_dm = descriptors[4];

                if (curName.equalsIgnoreCase("SpDiam") && (curMat.equals("EA")))
                    SpDiam_EA_dm = descriptors[5];
            }
        }
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
                int nCl = 0;
                for (IAtom at : CurMol.getConnectedAtomsList(at1)) {
                    if (at.getSymbol().equalsIgnoreCase("Cl"))
                        nCl++;
                }
                if (nCl == 1)
                    return 1.56;
                if (nCl == 2)
                    return 1.20;
                if (nCl == 3)
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
            int nH = 0;
            try {
                nH = at2.getImplicitHydrogenCount();
            } catch (Exception e) {
            }
            int nConn = CurMol.getConnectedBondsCount(at2) + nH;
            if ((ord == IBond.Order.SINGLE) && (nConn == 2))
                return 3.2;
            if ((ord == IBond.Order.SINGLE) && (nConn == 1))
                return 2.0;
//                return 0.3;
            if ((ord == IBond.Order.DOUBLE) && (nConn == 1))
                return 2.0;
        }


        // S-[O-]
        if ((a.equalsIgnoreCase("S")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nConn = CurMol.getConnectedBondsCount(at2);
            if ((ord == IBond.Order.SINGLE) && (nConn == 2))
                return 2.9;
        }


        // C(*)(*)-C(*)(*)(*) , C(*)(*)-C , CC(*)(*)(*)
        if ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("C"))) {

            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;

            int nH1 = 0, nH2 = 0;
            try {
                nH1 = at1.getImplicitHydrogenCount();
            } catch (Exception E) {
            }
            try {
                nH2 = at2.getImplicitHydrogenCount();
            } catch (Exception E) {
            }

            int nConn1 = CurMol.getConnectedBondsCount(at1) + nH1;
            int nConn2 = CurMol.getConnectedBondsCount(at2) + nH2;

            if ((nConn1 == 3) && (nConn2 == 4))
                return 0.68;
            if ((nConn1 == 3) && (nConn2 == 2))
                return 1.15;
            if ((nConn1 == 2) && (nConn2 == 4))
                return 1.48;
        }

        return 0;
    }

    private class EigenvalueBasedDescriptors {
        public HashMap<String, Double> Descriptors;

        public double[] Calculate(double[][] EigMat, int nBO) {

            Descriptors = new HashMap<>();

            // Calculates eigenvalues
            Matrix DataMatrix = new Matrix(EigMat);
            double[] eigenvalues;

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("unable_eigenvalue");
                double[] results = new double[5];
                Arrays.fill(results, MISSING_VALUE);
                return results;
            }

            // check on precision
            for (int i = 0; i < eigenvalues.length; i++) {
                if (Math.abs(eigenvalues[i]) < 0.000001)
                    eigenvalues[i] = 0;
            }

            // Eigenvalue based descriptors
            double EigMax = eigenvalues[0];
            double EigMin = eigenvalues[0];

            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
                if (val < EigMin)
                    EigMin = val;
            }
            double NormEigMax = EigMax / (double) nBO;
            double EigDiam = EigMax - EigMin;
            // Eigenvalue spectral moments
            double[] SpecMoments = new double[16];

            for (int i = 0; i < 16; i++) {
                SpecMoments[i] = 0;
                for (double val : eigenvalues) {
                    if (Math.abs(val) > 0)
                        SpecMoments[i] += Math.pow(val, (i + 1));
                }
                SpecMoments[i] = Math.log(1 + SpecMoments[i]);
            }

            double[] results = new double[6];
            results[0] = NormEigMax;
            results[1] = SpecMoments[14];
            results[2] = SpecMoments[2];
            results[3] = (eigenvalues.length - 7 >= 0) ? eigenvalues[eigenvalues.length - 7] : 0;
            results[4] = (eigenvalues.length - 5 >= 0) ? eigenvalues[eigenvalues.length - 5] : 0;
            results[5] = EigDiam;

            return results;

        }
    }

    public ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        String TYPE_D = "D";
        String TYPE_A = "A";
        String TYPE_P = "P";
        String TYPE_N = "N";
        String TYPE_L = "L";
        String TYPE_CYC = "Cyc";

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
                log.warn("unable_count_h");
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
            if (tA) AtomTypes[i].add(TYPE_A);
            if (tN) AtomTypes[i].add(TYPE_N);
            if (tP) AtomTypes[i].add(TYPE_P);
            if (tD) AtomTypes[i].add(TYPE_D);
            if (tL) AtomTypes[i].add(TYPE_L);
            if (tCyc) AtomTypes[i].add(TYPE_CYC);

        }

        return AtomTypes;
    }

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }

    private boolean isInCatsType(ArrayList<String> AtomA, ArrayList<String> AtomB, String TypeA, String TypeB) {
        if ( (isIn(AtomA, TypeA)) && (isIn(AtomB, TypeB)) )
            return true;
        if ( (isIn(AtomA, TypeB)) && (isIn(AtomB, TypeA)) )
            return true;
        return false;
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

    private int CalculateBin(String curWeight, double value) {

        String W_MASS = (new WeightsMass()).getSymbol();
        String W_ELECTRONEGATIVITY = (new WeightsElectronegativity()).getSymbol();

        if (curWeight.equalsIgnoreCase(W_MASS)) {
            // Mass
            if (value < 1.0) return 1;
            if (value < 1.2) return 2;
            if (value < 1.6) return 3;
            if (value < 3.0) return 4;
            return 5;
        }

        if (curWeight.equalsIgnoreCase(W_ELECTRONEGATIVITY)) {
            // Sanderson electronegativity
            if (value < 1.0) return 1;
            if (value < 1.1) return 2;
            if (value < 1.2) return 3;
            if (value < 1.3) return 4;
            if (value < 1.4) return 5;
            return 6;
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

        Object[][] RefBondLengths = {
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

        for (int i=0; i<RefBondLengths.length; i++) {
            if (AtomCouple(at1, at2, (String)RefBondLengths[i][0], (String)RefBondLengths[i][1])) {
                len = (Double)RefBondLengths[i][2];
                break;
            }
        }
        return len;
    }

    private class GCAtomCentredFragments {

        public final String[][] ACF_NAMES = {
                { "0", "U-000", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "1", "C-001", "CH3R / CH4"},
                { "2", "C-002", "CH2R2"},
                { "3", "C-003", "CHR3"},
                { "4", "C-004", "CR4"},
                { "5", "C-005", "CH3X"},
                { "6", "C-006", "CH2RX"},
                { "7", "C-007", "CH2X2"},
                { "8", "C-008", "CHR2X"},
                { "9", "C-009", "CHRX2"},
                { "10", "C-010", "CHX3"},
                { "11", "C-011", "CR3X"},
                { "12", "C-012", "CR2X2"},
                { "13", "C-013", "CRX3"},
                { "14", "C-014", "CX4"},
                { "15", "C-015", "=CH2"},
                { "16", "C-016", "=CHR"},
                { "17", "C-017", "=CR2"},
                { "18", "C-018", "=CHX"},
                { "19", "C-019", "=CRX"},
                { "20", "C-020", "=CX2"},
                { "21", "C-021", "#CH"},
                { "22", "C-022", "#CR / R=C=R"},
                { "23", "C-023", "#CX"},
                { "24", "C-024", "R--CH--R"},
                { "25", "C-025", "R--CR--R"},
                { "26", "C-026", "R--CX--R"},
                { "27", "C-027", "R--CH--X"},
                { "28", "C-028", "R--CR--X"},
                { "29", "C-029", "R--CX--X"},
                { "30", "C-030", "X--CH--X"},
                { "31", "C-031", "X--CR--X"},
                { "32", "C-032", "X--CX--X"},
                { "33", "C-033", "R--CH..X"},
                { "34", "C-034", "R--CR..X"},
                { "35", "C-035", "R--CX..X"},
                { "36", "C-036", "Al-CH=X"},
                { "37", "C-037", "Ar-CH=X"},
                { "38", "C-038", "Al-C(=X)-Al"},
                { "39", "C-039", "Ar-C(=X)-R"},
                { "40", "C-040", "R-C(=X)-X / R-C#X / X=C=X"},
                { "41", "C-041", "X-C(=X)-X"},
                { "42", "C-042", "X--CH..X"},
                { "43", "C-043", "X--CR..X"},
                { "44", "C-044", "X--CX..X"},
                { "45", "U-045", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "46", "H-046", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_46")},
                { "47", "H-047", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_47")},
                { "48", "H-048", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_48")},
                { "49", "H-049", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_49")},
                { "50", "H-050", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_50")},
                { "51", "H-051", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_51")},
                { "52", "H-052", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_52")},
                { "53", "H-053", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_53")},
                { "54", "H-054", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_54")},
                { "55", "H-055", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_55")},
                { "56", "O-056", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_56")},
                { "57", "O-057", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_57")},
                { "58", "O-058", "=O"},
                { "59", "O-059", "Al-O-Al"},
                { "60", "O-060", "Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X"},
                { "61", "O-061", "O--"},
                { "62", "O-062", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_62")},
                { "63", "O-063", "R-O-O-R"},
                { "64", "Se-064", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_64")},
                { "65", "Se-065", "=Se"},
                { "66", "N-066", "Al-NH2"},
                { "67", "N-067", "Al2-NH"},
                { "68", "N-068", "Al3-N"},
                { "69", "N-069", "Ar-NH2 / X-NH2"},
                { "70", "N-070", "Ar-NH-Al"},
                { "71", "N-071", "Ar-NAl2"},
                { "72", "N-072", "RCO-N< / >N-X=X"},
                { "73", "N-073", "Ar2NH / Ar3N / Ar2N-Al / R..N..R"},
                { "74", "N-074", "R#N / R=N-"},
                { "75", "N-075", "R--N--R / R--N--X"},
                { "76", "N-076", "Ar-NO2 / R--N(--R)--O / RO-NO"},
                { "77", "N-077", "Al-NO2"},
                { "78", "N-078", "Ar-N=X / X-N=X"},
                { "79", "N-079", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_79")},
                { "80", "U-080", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "81", "F-081", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_81")},
                { "82", "F-082", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_82")},
                { "83", "F-083", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_83")},
                { "84", "F-084", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_84")},
                { "85", "F-085", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_85")},
                { "86", "Cl-086", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_86")},
                { "87", "Cl-087", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_87")},
                { "88", "Cl-088", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_88")},
                { "89", "Cl-089", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_89")},
                { "90", "Cl-090", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_90")},
                { "91", "Br-091", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_91")},
                { "92", "Br-092", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_92")},
                { "93", "Br-093", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_93")},
                { "94", "Br-094", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_94")},
                { "95", "Br-095", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_95")},
                { "96", "I-096", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_96")},
                { "97", "I-097", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_97")},
                { "98", "I-098", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_98")},
                { "99", "I-099", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_99")},
                { "100", "I-100", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_100")},
                { "101", "F-101", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_101")},
                { "102", "Cl-102", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_102")},
                { "103", "Br-103", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_103")},
                { "104", "I-104", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_104")},
                { "105", "U-105", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "106", "S-106", "R-SH"},
                { "107", "S-107", "R2S / RS-SR"},
                { "108", "S-108", "R=S"},
                { "109", "S-109", "R-SO-R"},
                { "110", "S-110", "R-SO2-R"},
                { "111", "Si-111", ">Si<"},
                { "112", "B-112", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_112")},
                { "113", "U-113", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "114", "U-114", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
                { "115", "P-115", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_115")},
                { "116", "P-116", "R3-P=X"},
                { "117", "P-117", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_117")},
                { "118", "P-118", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_118")},
                { "119", "P-119", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_119")},
                { "120", "P-120", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_120")}
        };


        private final IAtomContainer CurMol;
        private final int nSK;
        private final double[][] ConnAugMatrix;
        private final boolean[] AtomAromatic;
        private final boolean ExplicitHydrogen;
        private final HashMap<Integer, Integer> NotMappedFragCount;

        private final int[] FragAtomId;


        public GCAtomCentredFragments(IAtomContainer Mol, boolean HasExplicitHydrogen) {

            // Init all variables
            this.CurMol = Mol;
            this.ExplicitHydrogen = HasExplicitHydrogen;
            nSK = Mol.getAtomCount();
            NotMappedFragCount = new HashMap<>();
            ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
            AtomAromatic = new boolean[nSK];
            for (int i=0; i<nSK; i++)
                AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

            // Init array of ACF id, all to MV
            FragAtomId = new int[nSK];
            for (int i : FragAtomId)
                i = Descriptor.MISSING_VALUE;

            // Check all atoms and try to assign to an ACF
            for (int AtomIdx=0; AtomIdx<nSK; AtomIdx++)
                ProcessAtom(AtomIdx);

        }


        public int[] GetACF() {
            return this.FragAtomId;
        }

        /**
         * Provide an hashmap with all the fragments (id, count) that are not directly mapped in the molecule. This is
         * needed to retrieve H fragments when a H-depleted molecule is processed.
         *
         * @return
         */
        public HashMap<Integer, Integer> getNotMappedFragCount() { return this.NotMappedFragCount; }


        private void ProcessAtom(int AtomIndex) {

            IAtom CurAt = CurMol.getAtom(AtomIndex);

            // if H skip - it will be calculated from the attached atom
            if (ExplicitHydrogen)
                if (ConnAugMatrix[AtomIndex][AtomIndex] == 1)
                    return;

            // if halogen skip - it will be calculated from the attached atom
            if (IsAtomHalogen(CurAt.getAtomicNumber()))
                return;


            int VD=0, nH=0, Charge=0;
            int nSingle=0, nDouble=0, nTriple=0, nArom=0;

            int sX=0, dX=0, tX=0, aX=0, asX=0;
            int sR=0, dR=0, tR=0, aR=0;
            int sAr=0, dAr=0, aAr=0;
            int sAl=0, dAl=0, aAl=0;
            int sAlx=0, dAlx=0, aAlx=0;

            // note: asX is aromatic single bond like in pyrrole:
            // cn[H]c, csc, coc

            int C_OxiNumber=0, C_Hybridazion=0, C_CX=0, C_CM=0;

            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex)
                    continue;
                if ( (ConnAugMatrix[AtomIndex][j]>0) && (ConnAugMatrix[j][j] != 1) ) {
                    VD++;
                    int Z = (int)ConnAugMatrix[j][j];
                    double b = ConnAugMatrix[AtomIndex][j];

                    if (b==1) {
                        nSingle++;
                        if (IsAtomElectronegative(Z))
                            sX++;
                        if (Z==6)
                            sR++;
                        if (AtomAromatic[j])
                            sAr++;
                        else
                        if (Z==6) sAl++;
                        else sAlx++;
                    }
                    if (b==2) {
                        nDouble++;
                        if (IsAtomElectronegative(Z))
                            dX++;
                        if (Z==6)
                            dR++;
                        if (AtomAromatic[j])
                            dAr++;
                        else
                        if (Z==6) dAl++;
                        else dAlx++;
                    }
                    if (b==3) {
                        nTriple++;
                        if (IsAtomElectronegative(Z))
                            tX++;
                        if (Z==6)
                            tR++;
                    }
                    if (b==1.5) {
                        nArom++;
                        if (IsAtomElectronegative(Z)) {

                            // checks if is a pyrrole-like aromatic single bond
                            int elNegVD=0, elNegH=0, elNegCharge;
                            for (int k=0; k<nSK; k++)
                                if ((ConnAugMatrix[j][k]>0) && (ConnAugMatrix[k][k]!=1)) {
                                    if (j==k) continue;
                                    elNegVD++;
                                }
                            try {
                                if (ExplicitHydrogen) {
                                    for (IAtom connAt : CurMol.getConnectedAtomsList(CurMol.getAtom(j)))
                                        if (connAt.getAtomicNumber() == 1)
                                            elNegH++;
                                } else {
                                    elNegH = CurMol.getAtom(j).getImplicitHydrogenCount();
                                }
                            } catch (Exception e) {
                                elNegH = 0;
                            }
                            try {
                                elNegCharge = CurMol.getAtom(j).getFormalCharge();
                            } catch (Exception e) {
                                elNegCharge = 0;
                            }
                            elNegVD += elNegH - elNegCharge;

                            boolean IsPyrroleLikeArom=false;

                            if ((Z==7) && (elNegVD==3))
                                IsPyrroleLikeArom = true;
                            if ((Z==8) && (elNegVD==2))
                                IsPyrroleLikeArom = true;
                            if ((Z==16) && (elNegVD==2) )
                                IsPyrroleLikeArom = true;

                            if (IsPyrroleLikeArom)
                                asX++;
                            else
                                aX++;

                        }
                        if (Z==6)
                            aR++;
                        if (AtomAromatic[j])
                            aAr++;
                        else
                        if (Z==6) aAl++;
                        else aAlx++;
                    }
                }
            }

            // counts H
            if (ExplicitHydrogen) {
                for (int j=0; j<nSK; j++) {
                    if (j==AtomIndex) continue;
                    if (ConnAugMatrix[j][AtomIndex] == 1)
                        if (ConnAugMatrix[j][j] == 1)
                            nH++;
                }
            } else {
                try {
                    nH = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }
            }

            // formal charge
            try {
                Charge = CurAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }

            // If Carbon, calculates oxidation number and hybridization
            if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                int c_VD=0;
                for (int j=0; j<nSK; j++) {
                    if (j==AtomIndex) continue;
                    if ( (ConnAugMatrix[j][AtomIndex]>0) && (ConnAugMatrix[j][j] != 1) ) {
                        c_VD++;

                        if (ConnAugMatrix[j][j] == 6) {
                            // search for -C-X
                            for (int k=0; k<nSK; k++) {
                                if (k==j) continue;
                                if ((ConnAugMatrix[k][j]>0) && (IsAtomElectronegative((int)ConnAugMatrix[k][k])))
                                    C_CX++;
                            }
                        } else if (IsAtomMetal((int)ConnAugMatrix[j][j])) {
                            C_CM++;
                        }
                    }
                }

                C_OxiNumber += (sX) + (dX * 2) + (tX * 3);
                C_OxiNumber += (asX);
                if (aX > 1)
                    C_OxiNumber += (aX * 1.5);
                else
                    C_OxiNumber += (aX * 2);

                C_Hybridazion = c_VD + nH - 1;
                if (((int)C_OxiNumber)!=C_OxiNumber)
                    C_OxiNumber = (int)C_OxiNumber + 1;
            }


            //// Search for proper fragment


            //// Hydrogen fragments ////////////////////////////////////////////////

            if (nH > 0) {

                int H_type = 0;

                if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                    boolean IsAlphaCarbon = false;

                    // Checks for alpha carbon

                    if ((nSingle > 0) && (nDouble==0) && (nTriple==0) && (nArom==0)) {

                        for (int j=0; j<nSK; j++) {
                            if (j==AtomIndex)
                                continue;

                            // -C
                            if (ConnAugMatrix[AtomIndex][j] > 0) {

                                if (ConnAugMatrix[j][j]==1)
                                    continue;

                                if ((ConnAugMatrix[AtomIndex][j]==1) && (ConnAugMatrix[j][j]==6)) {

                                    int nCdX=0, nCtX=0, nCaX=0;

                                    for (int k=0; k<nSK; k++) {
                                        if (k==j)
                                            continue;
                                        int Z = (int) ConnAugMatrix[k][k];
                                        if ((ConnAugMatrix[j][k]>0) && (IsAtomElectronegative(Z))) {
                                            if (ConnAugMatrix[j][k]==2)
                                                nCdX++;
                                            if (ConnAugMatrix[j][k]==3)
                                                nCtX++;
                                            if (ConnAugMatrix[j][k]==1.5)
                                                nCaX++;
                                        }
                                    }

                                    if (((nCdX + nCtX + nCaX) == 1) || (nCaX == 2))
                                        IsAlphaCarbon = true;

                                } else {

                                    IsAlphaCarbon = false; break;

                                }
                            }
                        }
                    }

                    if (IsAlphaCarbon) {
                        H_type=51;

                    } else {

                        // C0sp3 (no X attached to next C)
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==0) && (C_CM==0))
                            H_type=46;

                        // C1sp3, C0sp2
                        if ( ((C_OxiNumber==1) && (C_Hybridazion==3) && (C_CM==0))||
                                ((C_OxiNumber==0) && (C_Hybridazion==2) && (C_CM==0)))
                            H_type=47;

                        // C2sp3, C1sp2, C0sp
                        if ( ((C_OxiNumber==2) && (C_Hybridazion==3)) ||
                                ((C_OxiNumber==1) && (C_Hybridazion==2)) ||
                                ((C_OxiNumber==0) && (C_Hybridazion==1)) )
                            H_type=48;

                        // C3sp3, C2sp2, C2sp2, C3sp
                        if ( ((C_OxiNumber==3) && (C_Hybridazion==3)) ||
                                ((C_OxiNumber==2) && (C_Hybridazion==2)) ||
                                ((C_OxiNumber==3) && (C_Hybridazion==2)) ||
                                ((C_OxiNumber==3) && (C_Hybridazion==1)) )
                            H_type=49;

                        // C0sp3 with 1 X atom attached to next C
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==1) && (C_CM==0))
                            H_type=52;

                        // C0sp3 with 2 X atom attached to next C
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==2))
                            H_type=53;

                        // C0sp3 with 3 X atom attached to next C
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==3))
                            H_type=54;

                        // C0sp3 with 4 X atom attached to next C
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX>=4))
                            H_type=55;

                    }

                } else {

                    // H to heteroatom
                    H_type=50;

                }

                if (ExplicitHydrogen) {

                    // Sets the id for H atoms as they are explicit
                    for (int idxH=0; idxH<nSK; idxH++) {
                        if (idxH == AtomIndex) continue;
                        if ((ConnAugMatrix[AtomIndex][idxH] == 1) && (ConnAugMatrix[idxH][idxH] == 1) )
                            FragAtomId[idxH] = H_type;
                    }

                } else {

                    // Implicit H: counts the number of H and store the value in an hashmap, as these atoms
                    // can not be mapped directly on the molecule
                    int CountValue = nH;
                    if (NotMappedFragCount.containsKey(H_type))
                        CountValue += NotMappedFragCount.get(H_type);
                    NotMappedFragCount.put(H_type, CountValue);
                }
            }



            //// Halogen fragments /////////////////////////////////////////////////

            // Search for halogens attached to current atom
            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex) continue;
                if (ConnAugMatrix[j][AtomIndex]>0) {

                    // Cl
                    if (ConnAugMatrix[j][j] == 17) {

                        if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                            // attached to C3sp3
                            if ((C_OxiNumber==3) && (C_Hybridazion==3))
                            { FragAtomId[j] = 88;}

                        }
                    }
                }
            }

            //// Atom: N ///////////////////////////////////////////////////////////

            if (ConnAugMatrix[AtomIndex][AtomIndex] == 7) {

                // Checks for particular aromatic form
                boolean AromPyridineLike=false, AromPyrroleLike=false;
                if (nArom >= 2) {
                    if ((VD + nH - Charge) == 2 )
                        AromPyridineLike = true;
                    if ((VD + nH - Charge) == 3 )
                        AromPyrroleLike = true;
                }

                // Checks for particular N fragments linked to C, O or X


                int nO=0, sRCO=0, sXdX=0, sXvd2dX=0;

                for (int j=0; j<nSK; j++) {
                    if (j==AtomIndex)
                        continue;
                    if (ConnAugMatrix[AtomIndex][j]>0) {
                        int Z = (int)ConnAugMatrix[j][j];

                        // O
                        if (Z==8)
                            nO++;

                        // C for -RCO
                        if ( (Z==6) && (ConnAugMatrix[AtomIndex][j] == 1) )  {
                            for (int k=0; k<nSK; k++) {
                                if ((k==j) || (k==AtomIndex)) continue;
                                if ((ConnAugMatrix[k][j] == 2) && ( (ConnAugMatrix[k][k]==8) || (ConnAugMatrix[k][k]==16) ) ) {
                                    sRCO++;
                                    break;
                                }
                            }
                        }

                        if ((IsAtomElectronegative(Z))) {
                            int x_VD=0, x_dX=0;
                            for (int k=0; k<nSK; k++) {
//                            if ((k==j) || (k==AtomIndex)) continue;
                                if (k==j) continue;
                                if ((ConnAugMatrix[k][j]>0) && (ConnAugMatrix[k][k] != 1)) {
                                    x_VD++;
                                    if (IsAtomElectronegative((int)ConnAugMatrix[k][k]))
                                        if (ConnAugMatrix[k][j] == 2)
                                            x_dX++;
                                }
                            }
                            if ( (x_dX==1))
                                sXdX++;
                            if ( (x_dX==1) && (x_VD==2))
                                sXvd2dX++;
                        }

                    }
                }

                // fragments with possible [N+](=O)[O-] or R-O-N=O
                if (((VD==3) && (nH==0) && (sAr==1) && (nO==2)) ||
//                    ((VD==3) && (nH==0) && (aR==2) && (aX==1)) ||
                        ((VD==3) && (nH==0) && (aR==2) && (nO==1) && (AromPyridineLike)) ||
                        ((VD==2) && (nH==0) && (nO==2)))
//                || ((VD==2) && (nH==0) && (nO==2)))
                { FragAtomId[AtomIndex] = 76; return; }

                if ((VD==3) && (nH==0) && (sAr==0) && (nO>=2))
                { FragAtomId[AtomIndex] = 77; return; }

                // R=N- is matched also if it is R=NH
                if (((VD==1) && (nH==0) && (tR==1)) ||
                        ((VD==2) && (nH==0) && (dR==1) && (nSingle==1)) ||
                        ((VD==1) && (nH==1) && (dR==1)))
                { FragAtomId[AtomIndex] = 74; return; }

                if     (((VD>=2) && (nH==0) && (sAr==1) && (dX==1)) ||
                        ((VD>=2) && (nH==0) && (sX==1) && (dX==1)) ||
                        ((VD>=2) && (nH==0) && (dX==1) && (dR==1)) ||
                        ((VD>=2) && (nH==0) && (dX==2)) && (Charge == 0))
                { FragAtomId[AtomIndex] = 78; return; }

                // N Charged +1
                if ((!AromPyrroleLike) && (!AromPyridineLike))
                    if ((Charge==1) && (nTriple == 0))
                    { FragAtomId[AtomIndex] = 79; return; }

                // fragment with particular groups

                if (!AromPyrroleLike)
                    if ((((VD+nH)==3) && (sRCO>0)) ||   // RCO-N<
                            (((VD+nH)==3) && (sXvd2dX>0)))     // >N-X=X
                    { FragAtomId[AtomIndex] = 72; return; }


                if ((VD==1) && (nH==2) && (sAr==0) && (sX==0) && (nSingle==1))
                { FragAtomId[AtomIndex] = 66; return; }

                if ((VD==2) && (nH==1) && (sAr==0) && (nSingle==2))
                { FragAtomId[AtomIndex] = 67; return; }

                if ((VD==3) && (nH==0) && (sAr==0) && (aAr==0) && (nSingle==3))
                { FragAtomId[AtomIndex] = 68; return; }

                if ((VD==1) && (nH==2) && ( (sAr==1) || (sX==1) ))
                { FragAtomId[AtomIndex] = 69; return; }

                if ((VD==2) && (nH==1) && (sAl + sAlx ==1) && (sAr==1))
                { FragAtomId[AtomIndex] = 70; return; }

                if ((VD==3) && (nH==0) && (sAl + sAlx ==2) && (sAr==1))
                { FragAtomId[AtomIndex] = 71; return; }

                if (((VD==2) && (nH==1) && (sAr==2)) ||
                        ((VD==3) && (nH==0) && (sAr==3)) ||
                        ((VD==3) && (nH==0) && (sAr==2) && (sAl + sAlx ==1)) ||
                        ((AromPyrroleLike)))
                { FragAtomId[AtomIndex] = 73; return; }


                if (((AromPyridineLike)) ) // ||
                { FragAtomId[AtomIndex] = 75; return; }

                return;
            }

            //// Halogen ions //////////////////////////////////////////////////////

            // Cl
            if ((ConnAugMatrix[AtomIndex][AtomIndex] == 17) && (VD==0))
            { FragAtomId[AtomIndex] = 102;}

        }


        /**
         * Return true if an atom is electronegative, i.e. one of the following:
         * O, N, S, P, B, Si, Se or Halogens (F, Cl, Br, I)
         *
         * @param AtomicNumber
         * @return
         */
        private boolean IsAtomElectronegative(int AtomicNumber) {

            // O, N, S, P, B, Si, Se, halogen

            if ((AtomicNumber==7)||(AtomicNumber==8)||(AtomicNumber==15)||
                    (AtomicNumber==16)||(AtomicNumber==34)||(AtomicNumber==9)||
                    (AtomicNumber==5)||(AtomicNumber==14)||
                    (AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53)) {
                return true;
            }
            return false;
        }
        private boolean IsAtomMetal(int AtomicNumber) {

            // only Sn, Pb, Hg, As, Se, Ge

            if ((AtomicNumber==32) || (AtomicNumber==33) || (AtomicNumber==50) ||
                    (AtomicNumber==80) || (AtomicNumber==82)) {
                return true;
            }
            return false;
        }

        /**
         * Return true if an atom is Halogen (F, Cl, Br, I)
         *
         * @param AtomicNumber
         * @return
         */
        private boolean IsAtomHalogen(int AtomicNumber) {
            if ((AtomicNumber==9)||(AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53))
                return true;
            return false;
        }

        public String[][] getACF_NAMES() {
            return ACF_NAMES;
        }
    }

}
