package insilico.vapour_pressure.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.alerts.builders.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.iWeight;
import insilico.core.descriptor.blocks.weights.other.*;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
//import lombok.Data;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.vapour_pressure.descriptors.utils.BaryszMatrixCorrect;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertBlock;
import insilico.core.alerts.AlertList;

import java.util.*;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;

    public double H_050 = MISSING_VALUE;
    public double F02_F_F = MISSING_VALUE;
    public double SM02_AEAdm = MISSING_VALUE;
    public double SpMax_AEAed = MISSING_VALUE;
    public double ATS2m = MISSING_VALUE;
    public double P_VSA_PPP_A = MISSING_VALUE;
    public double NPerc = MISSING_VALUE;
    public double piID = MISSING_VALUE;
    public double CATS2D_07_DD = MISSING_VALUE;
    public double CATS2D_01_DN = MISSING_VALUE;
    public double N_072 = MISSING_VALUE;
    public double nAmideE_ = MISSING_VALUE;
    public double X3v = MISSING_VALUE;
    public double TI1_L = MISSING_VALUE;
    public double SpMaxA_Dzm = MISSING_VALUE;
    public double BB_SA12 = MISSING_VALUE;
    public double nCONN = MISSING_VALUE;
    public double CATS2D_07_DL = MISSING_VALUE;
    public double B03_F_F = MISSING_VALUE;
    public double SaaNH = MISSING_VALUE;
    public double ED_3 = MISSING_VALUE;
    public double CATS2D_09_DA = MISSING_VALUE;
    public double LogP = MISSING_VALUE;
    public double nRNH2 = MISSING_VALUE;
    public double F07_C_S = MISSING_VALUE;
    public double F10_C_O = MISSING_VALUE;
    public double CATS2D_09_AL = MISSING_VALUE;
    public double BB_SA31c = MISSING_VALUE;
    public double SpMin5_Bhm = MISSING_VALUE;
    public double X1Av = MISSING_VALUE;
    public double N_068 = MISSING_VALUE;
    public double nRNR2 = MISSING_VALUE;
    public double F09_O_O = MISSING_VALUE;
    public double IC1 = MISSING_VALUE;
    public double B01_C_N = MISSING_VALUE;
    public double nN_N = MISSING_VALUE;
    public double F10_C_Cl = MISSING_VALUE;


//    public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
//        InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION=true;
//        InsilicoMolecule mol_norm = (InsilicoMolecule) mol.Clone();
//        IAtomContainer ac = InsilicoMoleculeNormalization.Normalize(mol_norm.GetStructure());
//        String SMI = SmilesMolecule.GenerateSmiles(ac);
//        mol_norm.SetSMILESAndStructure(SMI, ac);
//        InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION=false;
//        CalculateAllDescriptors(mol_norm);
//    }

    public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws Exception {

        CalculateLaplaceMatrix(mol);
        CalculateInformationContent(mol);
        CalculateED(mol);
        CalculateEStateIndices(mol);
        CalculateLogP(mol);
        CalculateFunctionalGroups(mol);
        CalculateBenigniBossa(mol);
        CalculateConnectivityIndices(mol);
        CalculateIonizableGroups(mol);
        CalculateWalkAndPath(mol);
        CalculateConstitutional(mol);
        CalculateCATS2D(mol);
        CalculatePVSA(mol);
        CalculateAutoCorrelationFunction(mol);
        CalculateBurdenEigenvalues(mol);
        CalculateMatrices2D(mol);
        CalculateEdgeAdjacencyIndices(mol);
        CalculateAtomCentredFragments(mol);
        CalculateAtomPairs(mol);

    }

    private void CalculateInformationContent(InsilicoMolecule mol){
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("invalid_structure");
            return;
        }

        // Gets matrices
        int[][] TopoDistMat;
        try {
            TopoDistMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = curMol.getAtomCount();
        int[] VertexDeg = WeightsVertexDegree.getWeights(curMol, false);
        double[] IStates = (new WeightsIState()).getWeights(curMol, false);

        // rGes
        double[] ETIS = new double[nSK];
        for (int i=0; i<nSK; i++) {
            ETIS[i] = IStates[i];
            if (ETIS[i] != Descriptor.MISSING_VALUE) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    ETIS[i] += (IStates[i] - IStates[j]) / Math.pow(TopoDistMat[i][j] + 1, 2);
                }
            }
        }

        double[] s0k_Equal = new double[nSK];
        for (int i=0; i<nSK; i++)
            s0k_Equal[i] = ETIS[i];

        for (int i=0; i<nSK; i++) {
            if (s0k_Equal[i] != Descriptor.MISSING_VALUE) {
                int k = 1;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (Math.abs(s0k_Equal[i] - s0k_Equal[j]) < 0.001) {
                        k++;
                        s0k_Equal[j] = Descriptor.MISSING_VALUE;
                    }
                }
                s0k_Equal[i] = Descriptor.MISSING_VALUE;
            }
        }

        int[] DistDeg = new int[nSK];
        for (int i=0; i<nSK; i++) {
            DistDeg[i] = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if (TopoDistMat[i][j] > 0)
                    DistDeg[i] += TopoDistMat[i][j];
            }
        }

        int[] Equal = new int[nSK];
        for (int i=0; i<nSK; i++)
            Equal[i] = DistDeg[i];

        for (int i=0; i<nSK; i++) {
            if (Equal[i] != 0) {
                int k = 1;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (Equal[i] == Equal[j]) {
                        k++;
                        Equal[j] = 0;
                    }
                }
                Equal[i] = 0;
            }
        }

        double IVDE = 0;

        int[] VerDegCount = new int[10];
        for (int i=0; i<10; i++) VerDegCount[i] = 0;
        for (int i=0; i<nSK; i++) {
            VerDegCount[VertexDeg[i]]++;
        }

        for (int i=0; i<10; i++) {
            if (VerDegCount[i]>0)
                IVDE = IVDE - ( ((double)VerDegCount[i]) / ((double)nSK)) * Math.log( ((double)VerDegCount[i]) / ((double)nSK) );
        }

        double IVDEM = 0;

        for (int i=0; i<10; i++) {
            if (VerDegCount[i]>0)
                IVDEM = IVDEM - ( ( (double)VerDegCount[i] / (double)nSK ) * Log(2, ( (double)VerDegCount[i] / (double)nSK ) ) );
        }

        IAtomContainer curMolH;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            curMolH = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        // Gets matrices
        double[][] ConnMatH;
        int[][] TopoDistMatH;
        try {
            ConnMatH= ConnectionAugMatrix.getMatrix(curMolH);
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSKH = curMolH.getAtomCount();
        int[] VertexDegH = WeightsVertexDegree.getWeights(curMolH, false);

        int CurLag = 1;
        // Create belonging class for each atom(vertex)
        ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSKH);
        for (int i=0; i<nSKH; i++) {
            IAtom atStart =  curMolH.getAtom(i);
            ArrayList<String> CurNeig = new ArrayList<>();
            for (int j=0; j<nSKH; j++) {
                if (i==j) continue;
                if (TopoDistMatH[i][j] == CurLag) {
                    IAtom atEnd =  curMolH.getAtom(j);
                    ShortestPaths shortestPaths = new ShortestPaths(curMolH, atStart);
                    List<IAtom> sp = Arrays.asList(shortestPaths.atomsTo(atEnd));
                    StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                    for (int k=0; k<(sp.size()-1); k++) {
                        int a = curMolH.indexOf(sp.get(k));
                        int b = curMolH.indexOf(sp.get(k + 1));
                        if (ConnMatH[a][b] == 1)
                            bufPath.append("s");
                        if (ConnMatH[a][b] == 2)
                            bufPath.append("d");
                        if (ConnMatH[a][b] == 3)
                            bufPath.append("t");
                        if (ConnMatH[a][b] == 1.5)
                            bufPath.append("a");
                        bufPath.append(sp.get(k + 1).getAtomicNumber());
                        bufPath.append("(").append(VertexDegH[curMolH.indexOf(sp.get(k + 1))]).append(")");
                    }
                    CurNeig.add(bufPath.toString());
                }
            }
            Collections.sort(CurNeig);
            NeigList.add(CurNeig);
        }

        // Calculates equivalence classes
        ArrayList<ArrayList<String>> G = new ArrayList<>();
        ArrayList<Integer> Gn = new ArrayList<>();
        for (int i=0; i<nSKH; i++) {
            ArrayList<String> CurNeig = NeigList.get(i);
            boolean foundMatch = false;
            for (int k=0; k<G.size(); k++) {
                if (CompareNeigVector(CurNeig, G.get(k))) {
                    foundMatch = true;
                    int buf = Gn.get(k);
                    Gn.set(k, (buf+1));
                    break;
                }
            }
            if (!foundMatch) {
                G.add(CurNeig);
                Gn.add(1);
            }
        }

        // Calculate IC and CIC indices

        double ic=0;
        for (int i=0; i<Gn.size(); i++)
            ic += ((double)Gn.get(i)/nSKH) * (Math.log((double)Gn.get(i)/nSKH));
        ic = (-1.00 / Math.log(2)) * ic;

        IC1 =  ic;

    }
    private void CalculateED(InsilicoMolecule mol) throws Exception {

        try {
            SmartsPattern SA = SmartsPattern.create("[$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)][O,o][$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)]", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            ED_3 = 0;

            if (SA.matches(mol.GetStructure()))
                    ED_3 = 1;

        } catch (Throwable e) {
            log.warn("Unable to calculate" + e.getMessage());
        }
    }
    private void CalculateEStateIndices(InsilicoMolecule mol){
        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("invalid_structure");
            return;
        }

        int[][] TopoMatrix;
        try {
            TopoMatrix =  mol.GetMatrixTopologicalDistance();
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();

        // Get I-States weights
        WeightsIState w_istate = new WeightsIState();
        double[] w_is = w_istate.getWeights(m, false);
        for (double val : w_is)
            if (val == Descriptor.MISSING_VALUE) {
                log.warn("unable_calculate_all_estates");
                return;
            }

        double MaxDN = Descriptor.MISSING_VALUE;
        double MaxDP = Descriptor.MISSING_VALUE;

        // Calculate E-States
        double[] w_es = new double[nSK];
        for (int at = 0; at<nSK; at++) {
            double sumDeltaI = 0;

            for (int j = 0; j < nSK; j++)
                if (at != j)
                    sumDeltaI += (w_is[at] - w_is[j]) / Math.pow((double) TopoMatrix[at][j] + 1.0, 2.0);

            w_es[at] = w_is[at] + sumDeltaI;

            // max variation
            if (MaxDN == Descriptor.MISSING_VALUE)
                MaxDN = sumDeltaI;
            if (MaxDP == Descriptor.MISSING_VALUE)
                MaxDP = sumDeltaI;
            if (sumDeltaI > MaxDP)
                MaxDP = sumDeltaI;
            if (sumDeltaI < MaxDN)
                MaxDN = sumDeltaI;
        }

        SaaNH = 0;

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("unable_calculate_h_count");
            }

            // formal charge
            int Charge;
            try {
                Charge = curAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }

            // Count bonds
            int nSng = 0, nDbl = 0, nTri = 0, nAr = 0;
            for (IBond b : m.getConnectedBondsList(curAt)) {
                if (b.getFlag(CDKConstants.ISAROMATIC)) {
                    nAr++;
                    continue;
                }
                if (b.getOrder() == IBond.Order.SINGLE) {
                    nSng++;
                }
                if (b.getOrder() == IBond.Order.DOUBLE) {
                    nDbl++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nTri++;
                }
            }


            if (curAt.getSymbol().equals("N")) {
                if ((nSng == 0) && (nDbl == 0) &&
                        (nTri == 0) && (nAr == 2) &&
                        (nH == 1) && (Charge == 0)) {
                    SaaNH += w_es[at];
                }
            }
        }
    }
    private void CalculateLogP(InsilicoMolecule mol) throws InitFailureException, InvalidMoleculeException, GenericFailureException {

        SAMeylanLogPFragments Frags = new SAMeylanLogPFragments();
        Frags.Calculate(mol);
        SAMeylanLogPAdditionalFragments FragsAdd = new SAMeylanLogPAdditionalFragments();
        FragsAdd.Calculate(mol);
        SAMeylanLogPCorrectionFragments FragsCorr = new SAMeylanLogPCorrectionFragments();
        FragsCorr.Calculate(mol);

        double Coefficient = Frags.GetCoefficient() + FragsAdd.GetCoefficient();
        double Correction = FragsCorr.GetCoefficient();

        // main eq
        double Meylanlogp = 0.2290 + Coefficient + Correction;

        // lower bound
        if (Meylanlogp < -5.0)
            Meylanlogp = -5.0;

        LogP =  Meylanlogp;
    }
    private void CalculateFunctionalGroups(InsilicoMolecule mol) throws Exception {
        Pattern Query;
        String[] SMARTS = new String[]{
                "C(=[O,S])([$([#7;D3;!+](*)(*)*),$([#7;D2;!+](*)*),$([#7;D1;!+])])[$([#7;D3;!+](*)(*)*),$([#7;D2;!+](*)*),$([#7;D1;!+])]",
                "[N;D1][$([C;A]);!$(C=[O,S])]",
                "[N;D3]([$([C;A]);!$(C=[O,S])])([$([C;A]);!$(C=[O,S])])[$([C;A]);!$(C=[O,S])]",
                "[$([N;D1]),$([N;D2][C,c])]=[$([N;D1]),$([N;D2][C,c])]"
        };

        int nSMARTS = 0;

        for (String CurSMARTS : SMARTS) {
            try {
                Query = SmartsPattern.create(CurSMARTS).setPrepare(false);
            } catch (Exception e) {
                log.warn("descriptors_parsing_smarts_error");
                Query = null;
            }

            if (Query == null) {
                nSMARTS++;
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
                nSMARTS++;
                log.warn("descriptors_fgquery_init_fail" + e.getMessage());
                return;
            }

            if (nSMARTS == 0)
                nCONN = nmatch;

            if (nSMARTS == 1)
                nRNH2 = nmatch;

            if (nSMARTS == 2)
                nRNR2 = nmatch;

            if (nSMARTS == 3)
                nN_N = nmatch;

            nSMARTS++;
        }
    }
    private void CalculateBenigniBossa(InsilicoMolecule mol) throws InitFailureException, InvalidMoleculeException, GenericFailureException {
        AlertBlock BB;
        try {
            BB = new SABenigniBossa();
        } catch (Exception e) {
            BB = null;
            log.warn("Unable to init BB alert block - " + e.getMessage());
        }

        BB_SA12 = 0;
        BB_SA31c = 0;
        AlertList res = BB.Calculate(mol);
        for (Alert a: res.getSAList()){
            if (a.getId().equals("010c")){
                BB_SA12 = 1;
            }
            if (a.getId().equals("010z")){
                BB_SA31c = 1;
            }
        }
    }
    private void CalculateConnectivityIndices(InsilicoMolecule mol){
        ConnectivityIndices connectivityIndices = new ConnectivityIndices();
        connectivityIndices.Calculate(mol);
        X1Av = connectivityIndices.X1Av;
        X3v = connectivityIndices.X3v;
    }
    private void CalculateIonizableGroups(InsilicoMolecule mol) throws Exception {
        Pattern Query;
        String IG_SMARTS = "[!c][C](=[O])[$([N]([Br,Cl,I,F,O,N,P,S]));!$([N;R]([C;R]=[O,S]))]";

        try {
            Query = SmartsPattern.create(IG_SMARTS).setPrepare(false);
        } catch (Exception e) {
            log.warn("descriptors_parsing_smarts_error");
            Query = null;
        }

        if (Query == null)
            throw new Exception("descriptors_query_init_fail");

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

        nAmideE_ = nmatch;
    }
    private void CalculateWalkAndPath(InsilicoMolecule mol){
        try {
            MoleculePaths paths = new MoleculePaths(mol);
            double piID = paths.IDpi;
            piID = Math.log(1 + piID);
            this.piID = piID;

        } catch (GenericFailureException e) {
            log.warn("Error while calculating paths, unable to calculate");
        }
    }
    private void CalculateConstitutional(InsilicoMolecule mol){

        try {
            Constitutional constitutional = new Constitutional();
            constitutional.Calculate(mol);
            NPerc = constitutional.GetByName("NPerc").getValue();
        } catch (DescriptorNotFoundException e) {
            log.warn("invalid_structure");
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

        ArrayList<String>[] CatsTypes = setCatsAtomType(curMol, ConnAugMatrix);

        CATS2D_01_DN = 0;
        CATS2D_07_DL = 0;
        CATS2D_07_DD = 0;
        CATS2D_09_DA = 0;
        CATS2D_09_AL = 0;

        for (int i = 0; i < nSK; i++) {
            for (int j = i; j < nSK; j++) {

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "N"))
                    if (TopoMat[i][j] == 1)
                        CATS2D_01_DN++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "L"))
                    if (TopoMat[i][j] == 7)
                        CATS2D_07_DL++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "D"))
                    if (TopoMat[i][j] == 7)
                        CATS2D_07_DD++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "A"))
                    if (TopoMat[i][j] == 9)
                        CATS2D_09_DA++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "A", "L"))
                    if (TopoMat[i][j] == 9)
                        CATS2D_09_AL++;

            }
        }
    }
    private void CalculatePVSA(InsilicoMolecule mol){

        IAtomContainer curMolNoH;
        double[][] ConnAugMatrixNoH;
        try {
            curMolNoH = mol.GetStructure();
            ConnAugMatrixNoH = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn("invalid_structure");
            return;
        }

        int nSKnoH = curMolNoH.getAtomCount();

        ArrayList<String>[] AtomTypesOnMolWithoutH = setCatsAtomType(curMolNoH, ConnAugMatrixNoH);

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
            curDistance = Math.max(diffRadius, IdealBondLen[b]);

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

        double PVSA = 0;

        for (int i = 0; i < nSKnoH; i++) {

            for (String atType : AtomTypesOnMolWithoutH[i]){
                if (atType.equalsIgnoreCase("A")) {
                    PVSA += VSA[i];
                    break;
                }
            }
        }

        P_VSA_PPP_A =  PVSA;

    }
    private void CalculateAutoCorrelationFunction(InsilicoMolecule mol){
        iWeight curWeight = new WeightsMass();
        try {

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

            double[] w = ((iBasicWeight) curWeight).getScaledWeights(m);

            boolean MissingWeight = false;
            for (int i=0; i<nSK; i++)
                if (w[i] == Descriptor.MISSING_VALUE)
                    MissingWeight = true;
            if (MissingWeight)
                return;

            int lag=2;
            double AC=0;

            for (int i=0; i<nSK; i++) {
                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        AC += w[i] * w[j];
                    }
            }

            AC /= 2.0;
            AC = Math.log(1 + AC);


            ATS2m = AC;

        } catch (Throwable e) {
            log.warn("Unable to calculate" +  e.getMessage());
        }
    }
    private void CalculateBurdenEigenvalues(InsilicoMolecule mol){
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("invalid_structure");
            return;
        }

        // Gets matrix
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();

        iWeight curWeight = new WeightsMass();
        double[] w = ((iBasicWeight) curWeight).getScaledWeights(m);

        // If one or more weights are not available, sets all to missing value
        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE) {
                MissingWeight = true;
                break;
            }
        if (MissingWeight)
            return;

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        double valL;
        if (4 > (eigenvalues.length-1)) {
            valL = 0;
        } else {
            if (eigenvalues[4] < 0)
                valL = Math.abs(eigenvalues[4]);
            else
                valL = 0;
        }

        SpMin5_Bhm = valL;

    }

    private void CalculateLaplaceMatrix(InsilicoMolecule mol) throws GenericFailureException {
        try {
            IAtomContainer curMol;
            try {
                curMol = mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn("invalid_structure");
                return;
            }

            // Adj matrix available for calculations

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = new double[nSK][nSK];
            try {

                int[][] LapMatrix = mol.GetMatrixLaplace();
                for (int i=0; i<nSK; i++)
                    for (int j=0; j<nSK; j++)
                        Mat[i][j] = LapMatrix[i][j];

            } catch (Exception e) {
                log.warn(e.getMessage());
                return;
            }

            Matrix DataMatrix = new Matrix(Mat);
            double[] eigenvalues;

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("unable_eigenvalue" + e.getMessage());
                return;
            }

            // Laplace
            double QW;
            double TI1;
            double Sum = 0;
            double LastEigen;
            int matrix_dim = eigenvalues.length;

            for (int i=eigenvalues.length-1; i>=0; i--) {
                double val = eigenvalues[i];
                if (Math.abs(val) < Math.pow(1.5, -45))
                    continue;
                if (val > 0) {
                    LastEigen = val;
                    Sum += 1.0 / LastEigen;
                }
            }
            QW = Sum * matrix_dim;
            if (nBo > 0)
                TI1 = QW * 2 * Math.log((double)nBo / (double)matrix_dim) / Math.log(10);
            else
                TI1 = Descriptor.MISSING_VALUE;

            TI1_L = TI1;

        } catch (Throwable e) {
            log.warn("unable_calculate_matrices_2d" + e.getMessage());
        }

    }

    private void CalculateMatrices2D(InsilicoMolecule mol){
        try {

            IAtomContainer curMol;
            try {
                curMol = mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn("invalid_structure");
                return;
            }

            int nSK = curMol.getAtomCount();

            double[][] Mat = new double[nSK][nSK];
            try {
                double[][][] BarMat = BaryszMatrixCorrect.getMatrix(mol.GetStructure());
                int BarLayer = 1;

                for (int i=0; i<nSK; i++)
                    for (int j=0; j<nSK; j++)
                        Mat[i][j] = BarMat[i][j][BarLayer];
            } catch (Exception e) {
                log.warn(e.getMessage());
                return;
            }

            Matrix DataMatrix = new Matrix(Mat);
            double[] eigenvalues;

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("unable_eigenvalue" + e.getMessage());
                return;
            }

            double EigMax = eigenvalues[0];
            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
            }

            SpMaxA_Dzm = EigMax / (double) nSK;

        } catch (Throwable e) {
            log.warn("unable_calculate_matrices_2d" + e.getMessage());
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

        SpMax_AEAed = 0;

        double[][] curDataMatrix = new double[nBO][nBO];

        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        // weights for the matrix
        double[] w = new double[nBO];

        //// Edge degree
        for (int i=0; i<nBO; i++) {
            double deg = 0;
            for (int j=0; j<nBO; j++) {
                if (i==j) continue;
                if (curDataMatrix[i][j] != 0)
                    deg++;
            }
            w[i] = deg;
        }

        // AEA
        for (int i=0; i<nBO; i++)
            curDataMatrix[i][i] = w[i];

        EigenvalueBasedDescriptors eigDesc = new EigenvalueBasedDescriptors();
        SpMax_AEAed = eigDesc.Calculate(curDataMatrix, nBO)[0];

        //// Dipole moment
        for (int i=0; i<nBO; i++) {
            IAtom a =  curMol.getBond(i).getAtom(0);
            IAtom b =  curMol.getBond(i).getAtom(1);
            double CurVal = GetDipoleMoment(curMol, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(curMol, b, a);
            w[i] = CurVal;
        }

        // AEA
        for (int i=0; i<nBO; i++)
            curDataMatrix[i][i] = w[i];

        eigDesc = new EigenvalueBasedDescriptors();
        SM02_AEAdm = eigDesc.Calculate(curDataMatrix, nBO)[5];

    }
    private void CalculateAtomPairs(InsilicoMolecule mol){

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

        F02_F_F = 0;
        F07_C_S = 0;
        F09_O_O = 0;
        F10_C_O = 0;
        F10_C_Cl = 0;
        B01_C_N = 0;
        B03_F_F = 0;

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("F")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("F")) {
                        if (TopoMat[i][j] == 2)
                            F02_F_F++;
                        if (TopoMat[i][j] == 3)
                            B03_F_F = 1;
                    }
                }
            }
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("C")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("S")) {
                        if (TopoMat[i][j] == 7)
                            F07_C_S++;
                    }
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("O")) {
                        if (TopoMat[i][j] == 10)
                            F10_C_O++;
                    }
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("Cl")) {
                        if (TopoMat[i][j] == 10)
                            F10_C_Cl++;

                    }
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("N")) {
                        if (TopoMat[i][j] == 1)
                            B01_C_N = 1;
                    }
                }
            }
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("O")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("O")) {
                        if (TopoMat[i][j] == 9)
                            F09_O_O++;
                    }
                }
            }
        }
        F02_F_F /= 2;
        F09_O_O /= 2;
    }
    private void CalculateAtomCentredFragments(InsilicoMolecule mol){

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

        N_072 = FragDescriptors[72];
        N_068 = FragDescriptors[68];
        H_050 = FragDescriptors[50];

    }

    private class GCAtomCentredFragments {

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
            for (int i = 0; i < nSK; i++)
                AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

            // Init array of ACF id, all to MV
            FragAtomId = new int[nSK];
            for (int i : FragAtomId)
                i = Descriptor.MISSING_VALUE;

            // Check all atoms and try to assign to an ACF
            for (int AtomIdx = 0; AtomIdx < nSK; AtomIdx++)
                ProcessAtom(AtomIdx);

        }


        public int[] GetACF() {
            return this.FragAtomId;
        }

        public HashMap<Integer, Integer> getNotMappedFragCount() {
            return this.NotMappedFragCount;
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


            // note: asX is aromatic single bond like in pyrrole:
            // cn[H]c, csc, coc

            int C_OxiNumber = 0, C_Hybridazion = 0, C_CX = 0, C_CM = 0;

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
                        if (AtomAromatic[j])
                            sAr++;
                    }
                    if (b==2) {
                        nDouble++;
                        if (IsAtomElectronegative(Z))
                            dX++;
                        if (Z==6)
                            dR++;
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

                int H_type = 0;

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

                    if (IsAlphaCarbon) {
                        H_type = 51;

                    } else {

                        // C0sp3 (no X attached to next C)
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX == 0) && (C_CM == 0))
                            H_type = 46;

                        // C1sp3, C0sp2
                        if (((C_OxiNumber == 1) && (C_Hybridazion == 3) && (C_CM == 0)) ||
                                ((C_OxiNumber == 0) && (C_Hybridazion == 2) && (C_CM == 0)))
                            H_type = 47;

                        // C2sp3, C1sp2, C0sp
                        if (((C_OxiNumber == 2) && (C_Hybridazion == 3)) ||
                                ((C_OxiNumber == 1) && (C_Hybridazion == 2)) ||
                                ((C_OxiNumber == 0) && (C_Hybridazion == 1)))
                            H_type = 48;

                        // C3sp3, C2sp2, C2sp2, C3sp
                        if (((C_OxiNumber == 3) && (C_Hybridazion == 3)) ||
                                ((C_OxiNumber == 2) && (C_Hybridazion == 2)) ||
                                ((C_OxiNumber == 3) && (C_Hybridazion == 2)) ||
                                ((C_OxiNumber == 3) && (C_Hybridazion == 1)))
                            H_type = 49;

                        // C0sp3 with 1 X atom attached to next C
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX == 1) && (C_CM == 0))
                            H_type = 52;

                        // C0sp3 with 2 X atom attached to next C
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX == 2))
                            H_type = 53;

                        // C0sp3 with 3 X atom attached to next C
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX == 3))
                            H_type = 54;

                        // C0sp3 with 4 X atom attached to next C
                        if ((C_OxiNumber == 0) && (C_Hybridazion == 3) && (C_CX >= 4))
                            H_type = 55;

                    }

                } else {

                    // H to heteroatom
                    H_type = 50;

                }

                if (ExplicitHydrogen) {

                    // Sets the id for H atoms as they are explicit
                    for (int idxH = 0; idxH < nSK; idxH++) {
                        if (idxH == AtomIndex) continue;
                        if ((ConnAugMatrix[AtomIndex][idxH] == 1) && (ConnAugMatrix[idxH][idxH] == 1))
                            FragAtomId[idxH] = H_type;
                    }

                } else {

                    int CountValue = nH;
                    if (NotMappedFragCount.containsKey(H_type))
                        CountValue += NotMappedFragCount.get(H_type);
                    NotMappedFragCount.put(H_type, CountValue);
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


                int nO = 0, sRCO = 0, sXvd2dX = 0;

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
                            if ((x_dX == 1) && (x_VD == 2))
                                sXvd2dX++;
                        }

                    }
                }

                // fragments with possible [N+](=O)[O-] or R-O-N=O
                if (((VD == 3) && (nH == 0) && (sAr == 1) && (nO == 2)) ||
                        ((VD == 3) && (nH == 0) && (aR == 2) && (nO == 1) && (AromPyridineLike)) ||
                        ((VD == 2) && (nH == 0) && (nO == 2)))

                {
                    FragAtomId[AtomIndex] = 76;
                    return;
                }

                if ((VD == 3) && (nH == 0) && (sAr == 0) && (nO >= 2)) {
                    FragAtomId[AtomIndex] = 77;
                    return;
                }

                // R=N- is matched also if it is R=NH
                if (((VD == 1) && (nH == 0) && (tR == 1)) ||
                        ((VD == 2) && (nH == 0) && (dR == 1) && (nSingle == 1)) ||
                        ((VD == 1) && (nH == 1) && (dR == 1))) {
                    FragAtomId[AtomIndex] = 74;
                    return;
                }

//            System.out.println("AtomIndex = " + AtomIndex + ", Charge = " + Charge + ", sAr = " + sAr + ", dX = " + dX + ", sX = " + sX + ", dR = " + dR);
                if (((VD >= 2) && (nH == 0) && (sAr == 1) && (dX == 1)) ||
                        ((VD >= 2) && (nH == 0) && (sX == 1) && (dX == 1)) ||
                        ((VD >= 2) && (nH == 0) && (dX == 1) && (dR == 1)) ||
                        ((VD >= 2) && (nH == 0) && (dX == 2)) && (Charge == 0)) {
                    FragAtomId[AtomIndex] = 78;
                    return;
                }

                // N Charged +1
                if ((!AromPyrroleLike) && (!AromPyridineLike))
                    if ((Charge == 1) && (nTriple == 0)) {
                        FragAtomId[AtomIndex] = 79;
                        return;
                    }

                // fragment with particular groups

                if (!AromPyrroleLike)
                    if ((((VD + nH) == 3) && (sRCO > 0)) ||   // RCO-N<
                            (((VD + nH) == 3) && (sXvd2dX > 0)))     // >N-X=X
                    {
                        FragAtomId[AtomIndex] = 72;
                        return;
                    }

                if ((VD==1) && (nH==2) && (sAr==0) && (sX==0) && (nSingle==1))
                { FragAtomId[AtomIndex] = 66; return; }

                if ((VD==2) && (nH==1) && (sAr==0) && (nSingle==2))
                { FragAtomId[AtomIndex] = 67; return; }

                if ((VD==3) && (nH==0) && (sAr==0) && (aAr==0) && (nSingle==3))
                { FragAtomId[AtomIndex] = 68; }

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
    private class EigenvalueBasedDescriptors {

        public double[] Calculate(double[][] EigMat, int nBO) {

            Matrix DataMatrix = new Matrix(EigMat);
            double[] eigenvalues;

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("unable_eigenvalue" + e.getMessage());
                return new double[]{MISSING_VALUE,MISSING_VALUE,MISSING_VALUE,MISSING_VALUE,MISSING_VALUE};
            }

            // check on precision
            for (int i=0; i<eigenvalues.length; i++) {
                if (Math.abs(eigenvalues[i]) <  0.000001)
                    eigenvalues[i] = 0;
            }

            // Eigenvalue based descriptors
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

            EigAve = EigAve / (double) nBO;
            double NormEigMax = EigMax / (double) nBO;
            double EigDiam = EigMax - EigMin;
            double EigDev = 0;
            for (double val : eigenvalues)
                EigDev += Math.abs(val - EigAve);
            double NormEigDev = EigDev / (double) nBO;

            // Eigenvalue spectral moments
            double[] SpecMoments = new double[16];

            for (int i=0; i<2; i++) {
                SpecMoments[i] = 0;
                for (double val : eigenvalues) {
                    if (Math.abs(val) > 0)
                        SpecMoments[i] += Math.pow(val, (i + 1));
                }
                SpecMoments[i] = Math.log(1 + SpecMoments[i]);
            }

            return new double[]{EigMax, NormEigMax, EigDiam, EigDev, NormEigDev, SpecMoments[1]};

        }
    }
    public class MoleculePaths {

        // Limit for the number of atoms
        private static final int MAX_PATH_LENGTH = 2000;
        private static final int MAX_PATH_LENGTH_FOR_WALKS = 10;

        private IAtomContainer m;
        public double IDpi;
        private double[] Vertex_Distance_Degree;
        private boolean[] Entered;
        private double[] TotPCMult;
        private double[][] AdjConnectionMatrix;
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
            for (int i = 0; i < AdjMat.length; i++)
                for (int j = 0; j < AdjMat[0].length; j++)
                    AdjMatDbl[i][j] = AdjMat[i][j];
            double[][] AdjxAdj = new double[AdjMat.length][AdjMat[0].length];
            for (int i = 0; i < AdjMat.length; i++)
                for (int j = 0; j < AdjMat[0].length; j++)
                    AdjxAdj[i][j] = AdjMat[i][j];

            int[][] TopoDistMatrix;
            try {
                TopoDistMatrix = TopoDistanceMatrix.getMatrix(m);
            } catch (Exception e) {
                log.warn(e.getMessage());
                throw new GenericFailureException(e.getMessage());
            }
            try {
                AdjConnectionMatrix = Mol.GetMatrixConnectionAugmented();
            } catch (Exception e) {
                log.warn(e.getMessage());
                throw new GenericFailureException(e.getMessage());
            }


            nSK = m.getAtomCount();

            double[][] AdjPow = new double[nSK][nSK];

            double[][] AWC = new double[nSK][MAX_PATH_LENGTH];

            int Pow_k = 1;
            while (Pow_k < MAX_PATH_LENGTH_FOR_WALKS) {

                for (int i = 0; i < nSK; i++) {
                    for (int j = 0; j < nSK; j++) {
                        double ww = 0;
                        for (int m = 0; m < nSK; m++)
                            ww += AdjMatDbl[i][m] * AdjxAdj[m][j];
                        AdjPow[i][j] = ww;
                    }
                }

                for (int i = 0; i < nSK; i++) {
                    int aw = 0;
                    for (int j = 0; j < nSK; j++)
                        aw += AdjPow[i][j];
                    AWC[i][Pow_k] = aw;
                }

                Pow_k++;

                for (int i = 0; i < nSK; i++)
                    for (int j = 0; j < nSK; j++)
                        AdjxAdj[i][j] = AdjPow[i][j];

            }


            //// Paths calculation

            int LenChi = 5;

            int lim = Math.min(nSK, MAX_PATH_LENGTH);
            TotPCMult = new double[lim];

            Vertex_Distance_Degree = new double[nSK];

            for (int i = 0; i < nSK; i++) {
                Vertex_Distance_Degree[i] = 0;
                for (int j = 0; j < nSK; j++)
                    Vertex_Distance_Degree[i] += (double) TopoDistMatrix[i][j];
            }

            Entered = new boolean[nSK];
            int jlim = LenChi;

            for (int i = 0; i < nSK; i++) {

                for (int j = 0; j < nSK; j++)
                    Entered[j] = false;
                Entered[i] = true;

                int PathLength = 0;
                // Cycle on all atoms connected to i-th
                for (int j = 0; j < nSK; j++) {
                    if (j == i) continue;
                    if (AdjConnectionMatrix[i][j] != 0) {
                        double Cur_Mult_Bond_Order = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(i), m.getAtom(j)));
                        NextPathVisit(j, PathLength, Cur_Mult_Bond_Order);
                    }
                }
            }

            for (int i = 0; i < lim; i++) {
                TotPCMult[i] /= 2.0;
            }

            IDpi = nSK;

            for (int i = 0; i < lim; i++) {
                IDpi = IDpi + TotPCMult[i];
            }
        }


        private void NextPathVisit(int Atom_Idx, int PathLength, double Mult_Bond_Order){

            Entered[Atom_Idx] = true;
            if (PathLength < TotPCMult.length) {

                TotPCMult[PathLength] += Mult_Bond_Order;

                for (int Next_Atom_Idx = 0; Next_Atom_Idx < nSK; Next_Atom_Idx++) {
                    if (Next_Atom_Idx == Atom_Idx) continue;
                    if (Entered[Next_Atom_Idx]) continue;
                    if (AdjConnectionMatrix[Atom_Idx][Next_Atom_Idx] != 0) {

                        double BO = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(Atom_Idx), m.getAtom(Next_Atom_Idx)));

                        double Cur_Mult_Bond_Order = BO * Mult_Bond_Order;
                        NextPathVisit(Next_Atom_Idx, PathLength + 1, Cur_Mult_Bond_Order);

                    }
                }

            }
            Entered[Atom_Idx] = false;

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
    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }

    public ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        String TYPE_D = "D";
        String TYPE_A = "A";
        String TYPE_P = "P";
        String TYPE_N = "N";
        String TYPE_L = "L";

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
            if (tA) AtomTypes[i].add(TYPE_A);
            if (tN) AtomTypes[i].add(TYPE_N);
            if (tP) AtomTypes[i].add(TYPE_P);
            if (tD) AtomTypes[i].add(TYPE_D);
            if (tL) AtomTypes[i].add(TYPE_L);
//            if (tCyc) AtomTypes[i].add(TYPE_CYC[0]);

        }

        return AtomTypes;
    }
    private boolean isInCatsType(ArrayList<String> AtomA, ArrayList<String> AtomB, String TypeA, String TypeB) {
        if ( (isIn(AtomA, TypeA)) && (isIn(AtomB, TypeB)) )
            return true;
        if ( (isIn(AtomA, TypeB)) && (isIn(AtomB, TypeA)) )
            return true;
        return false;
    }
    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }
    private class ConnectivityIndices {

        private final static long serialVersionUID = 1L;
        private final static int MAX_PATH = 3;

        private boolean[] Entered;
        private double[] TotPC;
        private double[][] AdjConnectionMatrix;
        private int nSK;
        private int[] VD;
        private double[] ValenceVD;
        private int[] Qnumbers;
        private double[] curDescXv;
        double X3v = MISSING_VALUE;
        double X1Av = MISSING_VALUE;
        public void Calculate(InsilicoMolecule mol) {

            try {

                IAtomContainer curMol;
                try {
                    curMol = mol.GetStructure();
                } catch (InvalidMoleculeException e) {
                    log.warn("invalid_structure");
                    return;
                }

                // Gets matrices
                try {
                    AdjConnectionMatrix = mol.GetMatrixConnectionAugmented();
                } catch (GenericFailureException e) {
                    log.warn(e.getMessage());
                    return;
                }

                nSK = curMol.getAtomCount();

                VD = WeightsVertexDegree.getWeights(curMol, false);
                ValenceVD = (new WeightsValenceVertexDegree()).getWeights(curMol);
                Qnumbers = (new WeightsQuantumNumber()).getWeights(curMol);

                // checks for missing weights
                for (int i = 0; i < Qnumbers.length; i++)
                    if (Qnumbers[i] == -999)
                        return;

                for (int i = 0; i < ValenceVD.length; i++)
                    if (ValenceVD[i] == -999)
                        return;

                TotPC = new double[MAX_PATH];
                curDescXv = new double[MAX_PATH];

                Entered = new boolean[nSK];

                for (int k = 0; k < MAX_PATH; k++) {
                    curDescXv[k] = 0;
                }


                for (int i = 0; i < nSK; i++) {

                    for (int k = 0; k < nSK; k++)
                        Entered[k] = false;
                    Entered[i] = true;

                    int PathLength = 0;
                    boolean Atom1_F = (AdjConnectionMatrix[i][i] == 9);
                    int VD1 = VD[i];
                    int VD1_noHF = VD1;
                    for (int k = 0; k < nSK; k++) {
                        if (i == k) continue;
                        if ((AdjConnectionMatrix[i][k] > 0) && (AdjConnectionMatrix[k][k] == 9))
                            VD1_noHF--;
                    }

                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (AdjConnectionMatrix[i][j] > 0) {

                            boolean Atom2_F = (AdjConnectionMatrix[j][j] == 9);

                            int VD2 = VD[j];
                            int VD2_noHF = VD2;
                            for (int k = 0; k < nSK; k++) {
                                if (j == k) continue;
                                if ((AdjConnectionMatrix[j][k] > 0) && (AdjConnectionMatrix[k][k] == 9))
                                    VD2_noHF--;
                            }

                            double Chi_Solvation_Parameter;
                            if ((Atom1_F) || (Atom2_F))
                                Chi_Solvation_Parameter = 0;
                            else
                                Chi_Solvation_Parameter = (double) (Qnumbers[i] * Qnumbers[j]) / Math.sqrt(VD1_noHF * VD2_noHF);

                            NextPathVisit(j, PathLength, VD1 * VD2, ValenceVD[i] * ValenceVD[j],
                                    Chi_Solvation_Parameter);
                        }
                    }
                }

                double Xv;

                Xv = curDescXv[0] / 2.;
                X1Av = 0;
                if (TotPC[0] > 0)
                    X1Av = Xv / (TotPC[0] / 2.);

                X3v = curDescXv[2] / 2.;

            } catch (Throwable e) {
                log.warn("Unable to calculate:" + e.getMessage());
            }

        }
        private void NextPathVisit(int Atom_Idx, int PathLength, double Edge_Connec, double Valence_Edge_Connec,
                                   double Solvation) {

            Entered[Atom_Idx] = true;
            if (PathLength < curDescXv.length) {

                TotPC[PathLength] += 1;

                curDescXv[PathLength] += (1.0 / Math.sqrt(Valence_Edge_Connec));

                for (int Next_Atom_Idx = 0; Next_Atom_Idx < nSK; Next_Atom_Idx++) {
                    if (Next_Atom_Idx == Atom_Idx) continue;
                    if (Entered[Next_Atom_Idx]) continue;
                    if (AdjConnectionMatrix[Atom_Idx][Next_Atom_Idx] != 0) {

                        int VD2 = VD[Next_Atom_Idx];
                        int VD2_noHF = VD2;
                        for (int k = 0; k < nSK; k++) {
                            if (Next_Atom_Idx == k) continue;
                            if ((AdjConnectionMatrix[Next_Atom_Idx][k] > 0) && (AdjConnectionMatrix[k][k] == 9))
                                VD2_noHF--;
                        }

                        double SolvationParameter = 0;
                        if (VD2_noHF != 0)
                            SolvationParameter = Solvation * ((double) Qnumbers[Next_Atom_Idx] / Math.sqrt(VD2_noHF));

                        NextPathVisit(Next_Atom_Idx, PathLength + 1, Edge_Connec * VD2,
                                Valence_Edge_Connec * ValenceVD[Next_Atom_Idx],
                                SolvationParameter);
                    }
                }
            }
            Entered[Atom_Idx] = false;
        }
    }
    private boolean CompareNeigVector(ArrayList<String> A, ArrayList<String> B) {

        // note: integer in the vectors should be already sorted

        if (A.size()!=B.size())
            return false;
        else {
            for (int i=0; i<A.size(); i++)
                if (!A.get(i).equalsIgnoreCase(B.get(i)))
                    return false;
        }

        return true;
    }
    private double Log(int base,double x) {
        return Math.log10(x)/Math.log10((double)base);
    }
}
