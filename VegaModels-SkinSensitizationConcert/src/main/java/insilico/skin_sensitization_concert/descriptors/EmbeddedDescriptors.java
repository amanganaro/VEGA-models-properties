package insilico.skin_sensitization_concert.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertBlock;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.iWeight;
import insilico.core.descriptor.blocks.weights.other.*;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.*;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);
    private int MISSING_VALUE = -999;
    public double Eig12_AEA_ri = MISSING_VALUE;
    public double nArOR = MISSING_VALUE;
    public double nRCONH2 = MISSING_VALUE;
    public double nROH = MISSING_VALUE;
    public double CATS2D_08_DA = MISSING_VALUE;
    public double MATS3e = MISSING_VALUE;
    public double N_075 = MISSING_VALUE;
    public double F05_C_O = MISSING_VALUE;
    public double H_054 = MISSING_VALUE;
    public double nR_Cs = MISSING_VALUE;
    public double P_VSA_LogP_8 = MISSING_VALUE;
    public double CATS2D_09_AP = MISSING_VALUE;
    public double SRW03 = MISSING_VALUE;
    public double nSH = MISSING_VALUE;
    public double F06_N_F = MISSING_VALUE;
    public double B06_C_O = MISSING_VALUE;
    public double BB_SA18 = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws Exception {
        CalculateBenigniBossaAlerts(mol);
        CalculateWalkAndPathDescriptors(mol);
        CalculatePVSADescriptors(mol);
        Calculate2DAtomPairs(mol);
        CalculateACF(mol);
        Calculate2DAutoCorrelations(mol);
        CalculateCATS2D(mol);
        CalculateFunctionalGroupCount(mol);
        CalculateEdgeAdjacencyIndices(mol);
    }

    private void CalculateBenigniBossaAlerts(InsilicoMolecule mol) throws InvalidMoleculeException, GenericFailureException {
        AlertBlock BB;
        try {
            BB = new SABenigniBossa();
        } catch (Exception e) {
            BB = null;
            log.warn("Unable to init BB alert block - " + e.getMessage());
        }

        BB_SA18 = 0;
        AlertList res = BB.Calculate(mol);
        for (Alert a: res.getSAList()){
            if (a.getId().equals("010i")){
                BB_SA18 = 1;
            }
        }
    }
    private void CalculateWalkAndPathDescriptors(InsilicoMolecule mol){
        try {

            MoleculePaths paths = new MoleculePaths(mol);
            SRW03 = Math.log(1 + paths.Self_Returning_Walk_Counts[2]);

        } catch (GenericFailureException e) {
            log.warn("Error while calculating paths, unable to calculate: Walk and Path descriptors");
        }
    }
    private void CalculatePVSADescriptors(InsilicoMolecule mol){

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

        double[] w = null;
        w = (new WeightsHydrophobicityGC()).getWeightsForFragmentId(ACF);
        P_VSA_LogP_8 = 0;

        for (int i = 0; i < nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE) continue;
            if (CalculateBin(w[i]) == 8)
                P_VSA_LogP_8 += VSA[i];
        }
    }
    private void Calculate2DAtomPairs(InsilicoMolecule mol){
        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }
        int nSK = m.getAtomCount();

        int[][] TopoMat;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn("Invalid structure");
            return;
        }

        F05_C_O = 0;
        F06_N_F = 0;
        B06_C_O = 0;

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("C")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("O")) {
                        if (TopoMat[i][j] == 5)
                            F05_C_O++;

                        if (TopoMat[i][j] == 6)
                            B06_C_O = 1;
                    }
                }
            }
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("N")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("F")) {
                        if (TopoMat[i][j] == 6)
                            F06_N_F++;
                    }
                }
            }
        }
    }
    private void CalculateACF(InsilicoMolecule mol){
        IAtomContainer CurMol;
        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        new GCAtomCentredFragments(CurMol, false);
    }
    private void Calculate2DAutoCorrelations(InsilicoMolecule mol){
        try {

            IAtomContainer m;
            try {
                IAtomContainer orig_m = mol.GetStructure();
                m = Manipulator.AddHydrogens(orig_m);
            } catch (InvalidMoleculeException | GenericFailureException e) {
                log.warn("invalid_structure");
                return;
            }

            int[][] TopoMatrix;
            try {
                TopoMatrix = TopoDistanceMatrix.getMatrix(m);
            } catch (Exception e) {
                log.warn(e.getMessage());
                return;
            }

            int nSK = m.getAtomCount();

            double[] w;

            iWeight curWeight = new WeightsElectronegativity();
            w = ((iBasicWeight) curWeight).getScaledWeights(m);


            // If one or more weights are not available, sets all to missing value
            boolean MissingWeight = false;
            for (int i=0; i<nSK; i++)
                if (w[i] == Descriptor.MISSING_VALUE)
                    MissingWeight = true;
            if (MissingWeight)
                return;

            // Calculates weights averages
            double wA = 0;
            for (int i=0; i<nSK; i++)
                wA += w[i];
            wA = wA / ((double) nSK);


            double MoranAC=0;
            double denom = 0, delta = 0;

            for (int i=0; i<nSK; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == 3) {
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

            MATS3e = MoranAC;



        } catch (Throwable e) {
            log.warn("Unable to calculate: AutoCorrelation descriptors - " + e.getMessage());
        }
    }
    private void CalculateCATS2D(InsilicoMolecule mol){
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

        CATS2D_08_DA = 0;
        CATS2D_09_AP = 0;

        for (int i = 0; i < nSK; i++) {
            for (int j = i; j < nSK; j++) {

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "A"))
                    if (TopoMat[i][j] == 8 )
                        CATS2D_08_DA++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "A", "P"))
                    if (TopoMat[i][j] == 9 )
                        CATS2D_09_AP++;

            }
        }

    }
    private void CalculateFunctionalGroupCount(InsilicoMolecule mol) throws Exception {
        Pattern Query;
        String[] SMARTS = new String[] {
                "[#6;!$(C=[O,S]);!$(C#N)]O[c]",
                "[$([C;D3](=O)([N;D1;!+])[C;A]),$([C;D2](=O)[N;D1;!+])]",
                "[O;D1;!-]A",
                "[$([C;D2]([#6])=[#6]),$([C;D3]([!#6])([#6])=[#6])]",
                "[S;D1][C;!$(C=*);!$(C#*)]"
        };

        for (int smartsIdx = 0; smartsIdx < SMARTS.length; smartsIdx++) {
            String CurSMARTS = SMARTS[smartsIdx];

            try {
                Query = SmartsPattern.create(CurSMARTS).setPrepare(false);
            } catch (Exception e) {
                log.warn("descriptors_parsing_smarts_error");
                Query = null;
            }

            if (Query == null) {
                throw new Exception("descriptors_query_init_fail");
            }

            int nmatch = 0;

            List<Mappings> mappings = new ArrayList<>();
            boolean status;
            try {
                status = Query.matches(mol.GetStructure());
                if (status) {
                    mappings.add(Query.matchAll(mol.GetStructure()));
                    nmatch = Query.matchAll(mol.GetStructure()).countUnique();
                }
            } catch (Exception e) {
                log.warn("descriptors_fgquery_init_fail" + e.getMessage());
                return;
            }

            if (smartsIdx == 0)
                nArOR = nmatch;
            if (smartsIdx == 1)
                nRCONH2 = nmatch;
            if (smartsIdx == 2)
                nROH = nmatch;
            if (smartsIdx == 3)
                nR_Cs = nmatch;
            if (smartsIdx == 4)
                nSH = nmatch;
        }
    }
    private void CalculateEdgeAdjacencyIndices(InsilicoMolecule mol){
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

        // Gets basic matrix
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        double[][] curDataMatrix = new double[nBO][nBO];

        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        double[] w = new double[nBO];

        for (int i=0; i<nBO; i++)
            w[i] = GetResonanceIntegral(curMol.getBond(i));

        for (int i=0; i<nBO; i++)
            curDataMatrix[i][i] = w[i];


        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(curDataMatrix);
        double[] eigenvalues;

        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            log.warn("unable_eigenvalue" + e.getMessage());
            return;
        }

        // check on precision
        for (int i=0; i<eigenvalues.length; i++) {
            if (Math.abs(eigenvalues[i]) < 0.000001)
                eigenvalues[i] = 0;
        }

        int idx = eigenvalues.length - 12;
        double eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig12_AEA_ri = eig;

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
    private ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        ArrayList[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false, tCyc=false;

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


            // Sets final types
            if (tA) AtomTypes[i].add("A");
            if (tN) AtomTypes[i].add("N");
            if (tP) AtomTypes[i].add("P");
            if (tD) AtomTypes[i].add("D");
            if (tL) AtomTypes[i].add("L");

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

//        return 0.00;
        return 1.0;
    }
    private int CalculateBin(double value) {

            if (value < -1.5) return 1;
            if (value < -0.5) return 2;
            if (value < -0.25) return 3;
            if (value < 0) return 4;
            if (value < 0.25) return 5;
            if (value < 0.52) return 6;
            if (value < 0.75) return 7;
            return 8;

    }
    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }
    private double GetRefBondLength(IAtom at1, IAtom at2) {

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
    private class GCAtomCentredFragments {
        private final IAtomContainer CurMol;
        private final int nSK;
        private final double[][] ConnAugMatrix;
        private final boolean[] AtomAromatic;
        private final boolean ExplicitHydrogen;

        public GCAtomCentredFragments(IAtomContainer Mol, boolean HasExplicitHydrogen) {

            // Init all variables
            this.CurMol = Mol;
            this.ExplicitHydrogen = HasExplicitHydrogen;
            nSK = Mol.getAtomCount();
            ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
            AtomAromatic = new boolean[nSK];
            for (int i = 0; i < nSK; i++)
                AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

            N_075 = 0;
            H_054 = 0;

            // Check all atoms and try to assign to an ACF
            for (int AtomIdx = 0; AtomIdx < nSK; AtomIdx++)
                ProcessAtom(AtomIdx);

        }

        private void ProcessAtom(int AtomIndex) {

            IAtom CurAt = CurMol.getAtom(AtomIndex);

            // if H skip - it will be calculated from the attached atom
            if (ExplicitHydrogen)
                if (ConnAugMatrix[AtomIndex][AtomIndex] == 1)
                    return;

            // if halogen skip - it will be calculated from the attached atom
            if (IsAtomHalogen(CurAt.getAtomicNumber()))
                return;


            int VD = 0, nH = 0, Charge = 0;
            int nSingle = 0, nDouble = 0, nTriple = 0, nArom = 0;

            int sX = 0, dX = 0, tX = 0, aX = 0, asX = 0;
            int dR = 0, tR = 0, aR = 0;
            int sAr = 0, aAr = 0;
            int sAl = 0;
            int sAlx = 0;

            int C_OxiNumber = 0, C_Hybridazion = 0, C_CX = 0, C_CM = 0;

            for (int j = 0; j < nSK; j++) {
                if (j == AtomIndex)
                    continue;
                if ((ConnAugMatrix[AtomIndex][j] > 0) && (ConnAugMatrix[j][j] != 1)) {
                    VD++;
                    int Z = (int) ConnAugMatrix[j][j];
                    double b = ConnAugMatrix[AtomIndex][j];

                    if (b == 1) {
                        nSingle++;
                        if (IsAtomElectronegative(Z))
                            sX++;
                        if (AtomAromatic[j])
                            sAr++;
                        else if (Z == 6) sAl++;
                        else sAlx++;
                    }
                    if (b == 2) {
                        nDouble++;
                        if (IsAtomElectronegative(Z))
                            dX++;
                        if (Z == 6)
                            dR++;
                    }
                    if (b == 3) {
                        nTriple++;
                        if (IsAtomElectronegative(Z))
                            tX++;
                        if (Z == 6)
                            tR++;
                    }
                    if (b == 1.5) {
                        nArom++;
                        if (IsAtomElectronegative(Z)) {

                            // checks if is a pyrrole-like aromatic single bond
                            int elNegVD = 0, elNegH = 0, elNegCharge;
                            for (int k = 0; k < nSK; k++)
                                if ((ConnAugMatrix[j][k] > 0) && (ConnAugMatrix[k][k] != 1)) {
                                    if (j == k) continue;
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

                            boolean IsPyrroleLikeArom = false;

                            if ((Z == 7) && (elNegVD == 3))
                                IsPyrroleLikeArom = true;
                            if ((Z == 8) && (elNegVD == 2))
                                IsPyrroleLikeArom = true;
                            if ((Z == 16) && (elNegVD == 2))
                                IsPyrroleLikeArom = true;

                            if (IsPyrroleLikeArom)
                                asX++;
                            else
                                aX++;

                        }
                        if (Z == 6)
                            aR++;
                        if (AtomAromatic[j])
                            aAr++;
                    }
                }
            }

            // counts H
            if (ExplicitHydrogen) {
                for (int j = 0; j < nSK; j++) {
                    if (j == AtomIndex) continue;
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

                int c_VD = 0;
                for (int j = 0; j < nSK; j++) {
                    if (j == AtomIndex) continue;
                    if ((ConnAugMatrix[j][AtomIndex] > 0) && (ConnAugMatrix[j][j] != 1)) {
                        c_VD++;
                        if (ConnAugMatrix[j][j] == 6) {
                            // search for -C-X
                            for (int k = 0; k < nSK; k++) {
                                if (k == j) continue;
                                if ((ConnAugMatrix[k][j] > 0) && (IsAtomElectronegative((int) ConnAugMatrix[k][k])))
                                    C_CX++;
                            }
                        } else if (IsAtomMetal((int) ConnAugMatrix[j][j])) {
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
                if (((int) C_OxiNumber) != C_OxiNumber)
                    C_OxiNumber = (int) C_OxiNumber + 1;
            }

            //// Hydrogen fragments ////////////////////////////////////////////////

            if (nH > 0) {

                if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                    boolean IsAlphaCarbon = false;

                    // Checks for alpha carbon

                    if ((nSingle > 0) && (nDouble == 0) && (nTriple == 0) && (nArom == 0)) {

                        for (int j = 0; j < nSK; j++) {
                            if (j == AtomIndex)
                                continue;

                            // -C
                            if (ConnAugMatrix[AtomIndex][j] > 0) {

                                if (ConnAugMatrix[j][j] == 1)
                                    continue;

                                if ((ConnAugMatrix[AtomIndex][j] == 1) && (ConnAugMatrix[j][j] == 6)) {

                                    int nCdX = 0, nCtX = 0, nCaX = 0;

                                    for (int k = 0; k < nSK; k++) {
                                        if (k == j)
                                            continue;
                                        int Z = (int) ConnAugMatrix[k][k];
                                        if ((ConnAugMatrix[j][k] > 0) && (IsAtomElectronegative(Z))) {
                                            if (ConnAugMatrix[j][k] == 2)
                                                nCdX++;
                                            if (ConnAugMatrix[j][k] == 3)
                                                nCtX++;
                                            if (ConnAugMatrix[j][k] == 1.5)
                                                nCaX++;
                                        }
                                    }

                                    if (((nCdX + nCtX + nCaX) == 1) || (nCaX == 2))
                                        IsAlphaCarbon = true;

                                } else {

                                    IsAlphaCarbon = false;
                                    break;

                                }
                            }
                        }
                    }

                    if (!IsAlphaCarbon) {

                        // C0sp3 with 3 X atom attached to next C
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX == 3))
                            H_054 += nH;
                    }
                }
            }

            if (ConnAugMatrix[AtomIndex][AtomIndex] == 7) {

                // Checks for particular aromatic form
                boolean AromPyridineLike = false, AromPyrroleLike = false;
                if (nArom >= 2) {
                    if ((VD + nH - Charge) == 2)
                        AromPyridineLike = true;
                    if ((VD + nH - Charge) == 3)
                        AromPyrroleLike = true;
                }

                int nO = 0, sRCO = 0, sXdX = 0, sXvd2dX = 0;

                for (int j = 0; j < nSK; j++) {
                    if (j == AtomIndex)
                        continue;
                    if (ConnAugMatrix[AtomIndex][j] > 0) {
                        int Z = (int) ConnAugMatrix[j][j];

                        // O
                        if (Z == 8)
                            nO++;

                        // C for -RCO
                        if ((Z == 6) && (ConnAugMatrix[AtomIndex][j] == 1)) {
                            for (int k = 0; k < nSK; k++) {
                                if ((k == j) || (k == AtomIndex)) continue;
                                if ((ConnAugMatrix[k][j] == 2) && ((ConnAugMatrix[k][k] == 8) || (ConnAugMatrix[k][k] == 16))) {
                                    sRCO++;
                                    break;
                                }
                            }
                        }

                        if ((IsAtomElectronegative(Z))) {
                            int x_VD = 0, x_dX = 0;
                            for (int k = 0; k < nSK; k++) {
                                if (k == j) continue;
                                if ((ConnAugMatrix[k][j] > 0) && (ConnAugMatrix[k][k] != 1)) {
                                    x_VD++;
                                    if (IsAtomElectronegative((int) ConnAugMatrix[k][k]))
                                        if (ConnAugMatrix[k][j] == 2)
                                            x_dX++;
                                }
                            }
                            if ((x_dX == 1))
                                sXdX++;
                            if ((x_dX == 1) && (x_VD == 2))
                                sXvd2dX++;
                        }

                    }
                }

                if (((VD == 3) && (nH == 0) && (sAr == 1) && (nO == 2)) ||
                        ((VD == 3) && (nH == 0) && (aR == 2) && (nO == 1) && (AromPyridineLike)) ||
                        ((VD == 2) && (nH == 0) && (nO == 2))) {
                    return;
                }

                if ((VD == 3) && (nH == 0) && (sAr == 0) && (nO >= 2)) {
                    return;
                }

                // R=N- is matched also if it is R=NH
                if (((VD == 1) && (nH == 0) && (tR == 1)) ||
                        ((VD == 2) && (nH == 0) && (dR == 1) && (nSingle == 1)) ||
                        ((VD == 1) && (nH == 1) && (dR == 1))) {
                    return;
                }

                if (((VD >= 2) && (nH == 0) && (sAr == 1) && (dX == 1)) ||
                        ((VD >= 2) && (nH == 0) && (sX == 1) && (dX == 1)) ||
                        ((VD >= 2) && (nH == 0) && (dX == 1) && (dR == 1)) ||
                        ((VD >= 2) && (nH == 0) && (dX == 2)) && (Charge == 0)) {
                    return;
                }

                if (!AromPyrroleLike)
                    if ((((VD + nH) == 3) && (sRCO > 0)) ||
                            (((VD + nH) == 3) && (sXvd2dX > 0)))
                    {
                        return;
                    }

                if ((VD == 1) && (nH == 2) && (sAr == 0) && (sX == 0) && (nSingle == 1)) {
                    return;
                }

                if ((VD == 2) && (nH == 1) && (sAr == 0) && (nSingle == 2)) {
                    return;
                }

                if ((VD == 3) && (nH == 0) && (sAr == 0) && (aAr == 0) && (nSingle == 3)) {
                    return;
                }

                if ((VD == 1) && (nH == 2) && ((sAr == 1) || (sX == 1))) {
                    return;
                }

                if ((VD == 2) && (nH == 1) && (sAl + sAlx == 1) && (sAr == 1)) {
                    return;
                }

                if ((VD == 3) && (nH == 0) && (sAl + sAlx == 2) && (sAr == 1)) {
                    return;
                }

                if (((VD == 2) && (nH == 1) && (sAr == 2)) ||
                        ((VD == 3) && (nH == 0) && (sAr == 3)) ||
                        ((VD == 3) && (nH == 0) && (sAr == 2) && (sAl + sAlx == 1)) ||
                        ((AromPyrroleLike))) {
                    return;
                }


                if (AromPyridineLike) {
                    N_075++;
                    return;
                }
            }
        }

        private boolean IsAtomElectronegative(int AtomicNumber) {

            // O, N, S, P, B, Si, Se, halogen

            if ((AtomicNumber == 7) || (AtomicNumber == 8) || (AtomicNumber == 15) ||
                    (AtomicNumber == 16) || (AtomicNumber == 34) || (AtomicNumber == 9) ||
                    (AtomicNumber == 5) || (AtomicNumber == 14) ||
                    (AtomicNumber == 17) || (AtomicNumber == 35) || (AtomicNumber == 53)) {
                return true;
            }
            return false;
        }

        private boolean IsAtomMetal(int AtomicNumber) {

            // only Sn, Pb, Hg, As, Se, Ge

            if ((AtomicNumber == 32) || (AtomicNumber == 33) || (AtomicNumber == 50) ||
                    (AtomicNumber == 80) || (AtomicNumber == 82)) {
                return true;
            }
            return false;
        }

        private boolean IsAtomHalogen(int AtomicNumber) {
            if ((AtomicNumber == 9) || (AtomicNumber == 17) || (AtomicNumber == 35) || (AtomicNumber == 53))
                return true;
            return false;
        }
    }
    private class MoleculePaths {

        private static final int MAX_PATH_LENGTH = 2000;
        private static final int MAX_PATH_LENGTH_FOR_WALKS = 5;

        private IAtomContainer m;
        public double[] Self_Returning_Walk_Counts;
        private int nSK;



        /**
         * Constructor. When the object is created it directly calculates all path indices.
         * The input molecular structure should be provided WITHOUT explicit hydrogens.
         *
         * @param Mol
         * @throws GenericFailureException
         */
        public MoleculePaths(InsilicoMolecule Mol) throws GenericFailureException {
            Calculate(Mol);
        }


        private void Calculate(InsilicoMolecule Mol) throws GenericFailureException {

            try {
                m = Mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn("Invalid structure, unable to calculate paths");
                throw new GenericFailureException("invalid structure");
            }

            int[][] AdjMat;
            try {
                AdjMat = Mol.GetMatrixAdjacency();
            } catch (GenericFailureException e) {
                log.warn(e.getMessage());
                throw new GenericFailureException(e.getMessage());
            }
            double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
            for (int i=0; i<AdjMat.length; i++)
                for (int j=0; j<AdjMat[0].length; j++)
                    AdjMatDbl[i][j] = AdjMat[i][j];
            double[][] AdjxAdj = new double[AdjMat.length][AdjMat[0].length];
            for (int i=0; i<AdjMat.length; i++)
                for (int j=0; j<AdjMat[0].length; j++)
                    AdjxAdj[i][j] = AdjMat[i][j];

            nSK = m.getAtomCount();

            Self_Returning_Walk_Counts = new double[MAX_PATH_LENGTH];


            double[][] AdjPow = new double[nSK][nSK];

            Self_Returning_Walk_Counts[0] = nSK;

            int Pow_k = 1;
            while (Pow_k < MAX_PATH_LENGTH_FOR_WALKS) {

                for (int i=0; i<nSK; i++) {
                    for (int j = 0; j < nSK; j++) {
                        double ww = 0;
                        for (int m = 0; m < nSK; m++)
                            ww += AdjMatDbl[i][m] * AdjxAdj[m][j];
                        AdjPow[i][j] = ww;
                    }
                }

                for (int i=0; i<nSK; i++) {
                    Self_Returning_Walk_Counts[Pow_k] = Self_Returning_Walk_Counts[Pow_k] + AdjPow[i][i];
                }

                Pow_k++;

                for (int i=0; i<nSK; i++)
                    for (int j = 0; j < nSK; j++)
                        AdjxAdj[i][j] = AdjPow[i][j];

            }
        }
    }
}
