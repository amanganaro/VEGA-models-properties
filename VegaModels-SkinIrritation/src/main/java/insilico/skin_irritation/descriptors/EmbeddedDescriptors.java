package insilico.skin_irritation.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertBlock;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SABenigniBossa;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.Rings;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.iWeight;
import insilico.core.descriptor.blocks.weights.other.*;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.skin_irritation.descriptors.utils.BaryszMatrixCorrect;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
//import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.*;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;

    public double SM1_Dzi = MISSING_VALUE;
    public double Eig09_EAri = MISSING_VALUE;
    public double NdsCH = MISSING_VALUE;
    public double X3Av = MISSING_VALUE;
    public double NssS = MISSING_VALUE;
    public double nRNR2 = MISSING_VALUE;
    public double CATS2D_01_DA = MISSING_VALUE;
    public double IC3 = MISSING_VALUE;
    public double BB_SA44 = MISSING_VALUE;
    public double nArCHO = MISSING_VALUE;
    public double nRCO = MISSING_VALUE;
    public double SM03_AEAdm = MISSING_VALUE;
    public double B03_N_O = MISSING_VALUE;
    public double B04_C_N = MISSING_VALUE;
    public double Eig05_EAdm = MISSING_VALUE;
    public double C_026 = MISSING_VALUE;
    public double nCrs = MISSING_VALUE;
    public double SpMaxA_D_Dt = MISSING_VALUE;
    public double H_052 = MISSING_VALUE;
    public double nN = MISSING_VALUE;
    public double CATS2D_09_AL = MISSING_VALUE;
    public double Eig12_EAdm = MISSING_VALUE;
    public double Eig15_EAdm = MISSING_VALUE;
    public double CATS2D_08_NL = MISSING_VALUE;
    public double SpMin6_Bhs = MISSING_VALUE;
    public double SpMin6_Bhi = MISSING_VALUE;
    public double piPC10 = MISSING_VALUE;
    public double B06_O_O = MISSING_VALUE;
    public double Eig09_AEAri = MISSING_VALUE;
    public double SpMAD_AEAdm = MISSING_VALUE;
    public double B05_O_Cl = MISSING_VALUE;
    public double nR5 = MISSING_VALUE;
    public double nArOH = MISSING_VALUE;
    public double B04_Cl_Cl = MISSING_VALUE;
    public double IC1 = MISSING_VALUE;
    public double MPC09 = MISSING_VALUE;
    public double B09_O_O = MISSING_VALUE;
    public double CATS2D_08_PL = MISSING_VALUE;
    public double F08_C_N = MISSING_VALUE;
    public double GATS5s = MISSING_VALUE;
    public double F03_C_Cl = MISSING_VALUE;
    public double H_051 = MISSING_VALUE;
    public double B06_C_N = MISSING_VALUE;
    public double GATS1e = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws Exception {


        CalculateWalkAndPath(mol);
        CalculateBurdenEigenvalues(mol);
        CalculateConstitutionalIndices(mol);
        CalculateACF(mol);
        CalculateAtomPairsDescriptors(mol);
        CalculateICDescriptors(mol);
        CalculateCATS2DDescriptors(mol);
        CalculateFunctionalGroupsDescriptors(mol);
        CalculateConnectivityIndices(mol);
        CalculateEstateIndices(mol);
        CalculateEdgeAdjacencyIndices(mol);
        CalculateMatrixBasedDescriptors(mol);
    }

    private void CalculateBenigniBossaAlerts(InsilicoMolecule mol) throws InvalidMoleculeException, GenericFailureException {
        AlertBlock BB;
        try {
            BB = new SABenigniBossa();
        } catch (Exception e) {
            BB = null;
            log.warn("Unable to init BB alert block - " + e.getMessage());
        }

        BB_SA44 = 0;
        AlertList res = BB.Calculate(mol);
        for (Alert a: res.getSAList()){
            if (a.getId().equals("")){
                BB_SA44 = 1;
            }
        }
    }
    private void CalculateWalkAndPath(InsilicoMolecule mol){
        try {

            MoleculePaths paths = new MoleculePaths(mol);

            MPC09 = Math.log(1 + paths.Path_Counts[8]);
            piPC10 = Math.log(1 + paths.Multiple_Path_Counts[9]);


        } catch (GenericFailureException e) {
            log.warn("Error while calculating paths, unable to calculate Walk and Path descriptors");
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

        ArrayList<iWeight> bWeights = new ArrayList<>();
        bWeights.add(new WeightsIonizationPotential());
        bWeights.add(new WeightsIState());

        // Gets matrix
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();

        for (int weightIdx = 0; weightIdx<bWeights.size(); weightIdx++){

            iWeight curWeight = bWeights.get(weightIdx);

            double[] w;

            if (curWeight.getClass() == WeightsIState.class) {
                w = ((WeightsIState)curWeight).getWeights(m, true);
                for (int i=0; i<nSK; i++) {
                    if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                        w[i] = 1;
                }
            } else {
                w = ((iBasicWeight) curWeight).getScaledWeights(m);
            }

            boolean MissingWeight = false;
            for (int i=0; i<nSK; i++)
                if (w[i] == Descriptor.MISSING_VALUE) {
                    MissingWeight = true;
                    break;
                }
            if (MissingWeight)
                continue;

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
                if (eigenvalues.length < 6) {
                    valL = 0;
                } else {
                    if (eigenvalues[5] < 0)
                        valL = Math.abs(eigenvalues[5]);
                    else
                        valL = 0;
                }

                if (weightIdx == 0)
                    SpMin6_Bhi = valL;
                if (weightIdx == 1)
                    SpMin6_Bhs = valL;
        }
    }
    private void CalculateConstitutionalIndices(InsilicoMolecule mol) throws DescriptorNotFoundException {
        Constitutional constitutional = new Constitutional();
        constitutional.Calculate(mol);
        nN = constitutional.GetByName("nN").getValue();
    }
    private void CalculateACF(InsilicoMolecule mol){
        IAtomContainer CurMol;
        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return;
        }

        GCAtomCentredFragments GC = new GCAtomCentredFragments(CurMol, false);
        C_026 = GC.C_026;
        H_051 = GC.H_051;
        H_052 = GC.H_052;
    }
    private void CalculateAtomPairsDescriptors(InsilicoMolecule mol){

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

        B03_N_O = 0;
        B04_C_N = 0;
        B04_Cl_Cl = 0;
        B05_O_Cl = 0;
        B06_C_N = 0;
        B06_O_O = 0;
        B09_O_O = 0;
        F03_C_Cl = 0;
        F08_C_N = 0;

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("N")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("O")) {
                        if (TopoMat[i][j] == 3)
                            B03_N_O = 1;
                    }
                }
            }

            if (m.getAtom(i).getSymbol().equalsIgnoreCase("C")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("Cl")) {
                        if (TopoMat[i][j] == 3)
                            F03_C_Cl++;

                    }
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("N")) {
                        if (TopoMat[i][j] == 4) {
                            B04_C_N = 1;
                            B06_C_N = 1;
                        }
                        if (TopoMat[i][j] == 8)
                            F08_C_N++;
                    }
                }
            }

            if (m.getAtom(i).getSymbol().equalsIgnoreCase("O")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("O")) {
                        if (TopoMat[i][j] == 6)
                            B06_O_O = 1;
                        if (TopoMat[i][j] == 9)
                            B09_O_O = 1;
                    }
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("Cl")) {
                        if (TopoMat[i][j] == 5)
                            B05_O_Cl = 1;
                    }
                }
            }

            if (m.getAtom(i).getSymbol().equalsIgnoreCase("Cl")) {
                for (int j=0; j<nSK; j++) {
                    if (i == j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase("Cl")) {
                        if (TopoMat[i][j] == 4)
                            B04_Cl_Cl = 1;
                    }
                }
            }
        }
    }
    private void CalculateICDescriptors(InsilicoMolecule mol){
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

        for (int CurLag=1; CurLag<=3; CurLag+=2) {

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


            if (CurLag == 1)
                IC1 = ic;
            if (CurLag == 3)
                IC3 = ic;
        }
    }
    private void CalculateCATS2DDescriptors(InsilicoMolecule mol){
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

        CATS2D_01_DA = 0;
        CATS2D_08_PL = 0;
        CATS2D_08_NL = 0;
        CATS2D_09_AL = 0;

        int[] desc = new int[10];
        Arrays.fill(desc, 0);

        for (int i = 0; i < nSK; i++) {
            for (int j = i; j < nSK; j++) {

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "D", "A"))
                    if (TopoMat[i][j] == 1 )
                        CATS2D_01_DA++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "P", "L"))
                    if (TopoMat[i][j] == 8 )
                        CATS2D_08_PL++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "N", "L"))
                    if (TopoMat[i][j] == 8 )
                        CATS2D_08_NL++;

                if (isInCatsType(CatsTypes[i], CatsTypes[j], "A", "L"))
                    if (TopoMat[i][j] == 9 )
                        CATS2D_09_AL++;
            }
        }
    }
    private void CalculateFunctionalGroupsDescriptors(InsilicoMolecule mol) throws Exception {
        Pattern Query;

        double[] FuncGroupDescr = new double[]{nArOH, nCrs, nRCO, nArCHO, nRNR2};
        String[] SMARTS = new String[]{
                "[O;D1;!-]a",
                "[$([C;D2;H2](@[#6])@[#6]),$([C;D3;H](@[#6])(@[#6])[*;!#6]),$([C;D3;H](@[#6])(@[!#6])[#6]),$([C;D4](@[#6])(@[#6])([*;!#6])[*;!#6]),$([C;D4](@[#6])(@[!#6])([#6])[*;!#6]),$([C;D4](@[!#6])(@[!#6])([#6])[#6])]",
                "[C;D3](=O)([C])[C]",
                "[$([C;D2](=O)[c])]",
                "[N;D3]([$([C;A]);!$(C=[O,S])])([$([C;A]);!$(C=[O,S])])[$([C;A]);!$(C=[O,S])]"
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
                nArOH = nmatch;
            if (smartsIdx == 1)
                nCrs = nmatch;
            if (smartsIdx == 2)
                nRCO = nmatch;
            if (smartsIdx == 3)
                nArCHO = nmatch;
            if (smartsIdx == 4)
                nRNR2 = nmatch;

        }
    }
    private void CalculateConnectivityIndices(InsilicoMolecule mol){
        ConnectivityIndices CI = new ConnectivityIndices();
        CI.Calculate(mol);
        X3Av = CI.X3Av;
    }
    private void CalculateEstateIndices(InsilicoMolecule mol){
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
        for (int at = 0; at<nSK; at++) {
            double sumDeltaI = 0;

            for (int j = 0; j < nSK; j++)
                if (at != j)
                    sumDeltaI += (w_is[at] - w_is[j]) / Math.pow((double) TopoMatrix[at][j] + 1.0, 2.0);


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

        NdsCH = 0;
        NssS = 0;

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
            int nSng = 0, nDbl = 0, nTri = 0, nAr=0;
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

            if  (curAt.getSymbol().equals("C")) {
                if ( (nSng == 1)  && (nDbl == 1) &&
                        (nTri == 0) && (nAr == 0) &&
                        (nH == 1) && (Charge == 0)) {
                    NdsCH++;
                }
            }

            if  (curAt.getSymbol().equals("S")) {
                if ( (nSng == 2)  && (nDbl == 0) &&
                        (nTri == 0) && (nAr == 0) &&
                        (nH == 0) && (Charge == 0)) {
                    NssS++;
                }
            }
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

        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        // weights for the matrix
        double[] w = new double[nBO];

        /////////////////////// PARTE COMUNE!!

        //// Dipole moment

        for (int i=0; i<nBO; i++) {
            IAtom a =  curMol.getBond(i).getAtom(0);
            IAtom b =  curMol.getBond(i).getAtom(1);
            double CurVal = GetDipoleMoment(curMol, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(curMol, b, a);
            w[i] = CurVal;
        }

        for (int i=0; i<nBO; i++)
            curDataMatrix[i][i] = w[i];

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

        for (int i=0; i<eigenvalues.length; i++) {
            if (Math.abs(eigenvalues[i]) < 0.000001)
                eigenvalues[i] = 0;
        }

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
        double EigDev = 0;
        for (double val : eigenvalues)
            EigDev += Math.abs(val - EigAve);
        double NormEigDev = EigDev / (double) nBO;

        SpMAD_AEAdm =  NormEigDev;
        SM03_AEAdm = 0;
        for (double val : eigenvalues) {
            if (Math.abs(val) > 0)
                SM03_AEAdm += Math.pow(val, 3);
        }
        SM03_AEAdm = Math.log(1 + SM03_AEAdm);


        ///////////////////////////////

        curDataMatrix = new double[nBO][nBO];

        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<nBO; i++) {
            IAtom a =  curMol.getBond(i).getAtom(0);
            IAtom b =  curMol.getBond(i).getAtom(1);
            double CurVal = GetDipoleMoment(curMol, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(curMol, b, a);
            w[i] = CurVal;
        }

        for (int i=0; i<nBO; i++) {
            for (int j=0; j<nBO; j++) {
                if (curDataMatrix[i][j] != 0)
                    curDataMatrix[i][j] = w[j];
            }
        }

        DataMatrix = new Matrix(curDataMatrix);

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

        EigAve = 0;
        EigMax = eigenvalues[0];
        EigMin = eigenvalues[0];

        for (double val : eigenvalues) {
            EigAve += val;
            if (val > EigMax)
                EigMax = val;
            if (val< EigMin)
                EigMin = val;
        }
        EigAve = EigAve / (double) nBO;

        EigDev = 0;
        for (double val : eigenvalues)
            EigDev += Math.abs(val - EigAve);

        int idx = eigenvalues.length - 5;
        double eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig05_EAdm =  eig;

        idx = eigenvalues.length - 12;
        eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig12_EAdm =  eig;

        idx = eigenvalues.length - 15;
        eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig15_EAdm =  eig;


        ///////////////////////////////

        curDataMatrix = new double[nBO][nBO];

        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<nBO; i++)
            w[i] = GetResonanceIntegral(curMol.getBond(i));


        for (int i=0; i<nBO; i++) {
            for (int j=0; j<nBO; j++) {
                if (curDataMatrix[i][j] != 0)
                    curDataMatrix[i][j] = w[j];
            }
        }

        DataMatrix = new Matrix(curDataMatrix);

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

        EigAve = 0;
        EigMax = eigenvalues[0];
        EigMin = eigenvalues[0];

        for (double val : eigenvalues) {
            EigAve += val;
            if (val > EigMax)
                EigMax = val;
            if (val< EigMin)
                EigMin = val;
        }
        EigAve = EigAve / (double) nBO;

        EigDev = 0;
        for (double val : eigenvalues)
            EigDev += Math.abs(val - EigAve);


        idx = eigenvalues.length - 9;
        eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig09_EAri =  eig;

        ///////////////////////////////

        curDataMatrix = new double[nBO][nBO];

        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<nBO; i++)
            w[i] = GetResonanceIntegral(curMol.getBond(i));


        for (int i=0; i<nBO; i++)
            curDataMatrix[i][i] = w[i];

        DataMatrix = new Matrix(curDataMatrix);

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

        EigAve = 0;
        EigMax = eigenvalues[0];
        EigMin = eigenvalues[0];

        for (double val : eigenvalues) {
            EigAve += val;
            if (val > EigMax)
                EigMax = val;
            if (val< EigMin)
                EigMin = val;
        }
        EigAve = EigAve / (double) nBO;

        EigDev = 0;
        for (double val : eigenvalues)
            EigDev += Math.abs(val - EigAve);


        idx = eigenvalues.length - 9;
        eig = (idx >= 0) ? eigenvalues[idx] : 0;
        Eig09_AEAri =  eig;

    }
    private void CalculateMatrixBasedDescriptors(InsilicoMolecule mol){
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("invalid_structure");
            return;
        }

        // Adj matrix available for calculations
        int nSK = curMol.getAtomCount();

        String[] MATRICES = {"D/Dt", "Dz(i)"};

        for (String curMat : MATRICES) {

            // Gets current matrix
            double[][] Mat = new double[nSK][nSK];
            try {

                if (curMat.equalsIgnoreCase("D/Dt")) {
                    Mat = mol.GetMatrixDistanceDetour();
                }

                if (curMat.equalsIgnoreCase("Dz(i)")) {
                    double[][][] BarMat = BaryszMatrixCorrect.getMatrix(mol.GetStructure());
                    int BarLayer = 5;
                    for (int i = 0; i < nSK; i++)
                        for (int j = 0; j < nSK; j++)
                            Mat[i][j] = BarMat[i][j][BarLayer];
                }

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
            double EigMin = eigenvalues[0];
            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
                if (val< EigMin)
                    EigMin = val;
            }

            double NormEigMax = EigMax / (double) nSK;

            if (curMat.equals("Dz(i)")) {
                SM1_Dzi = 0;

                for (double val : eigenvalues)
                    SM1_Dzi += val;

                SM1_Dzi = Math.signum(SM1_Dzi) * Math.log(1 + Math.abs(SM1_Dzi));
            }

            if (curMat.equals("D/Dt"))
                SpMaxA_D_Dt =  NormEigMax;

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
    private class ConnectivityIndices {

        private final static int MAX_PATH = 3;
        private boolean[] Entered;
        private double[] TotPC;
        private double[][] AdjConnectionMatrix;
        private int nSK;
        private int[] VD;
        private double[] ValenceVD;
        private int[] Qnumbers;
        private double[] curDescXv;
        double X3Av = MISSING_VALUE;
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

                double Xv = curDescXv[2] / 2.;
                X3Av = 0;
                if (TotPC[2] > 0)
                    X3Av = Xv / (TotPC[2] / 2.);

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
    private ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        ArrayList[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false;

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
    private class GCAtomCentredFragments {

        private final IAtomContainer CurMol;
        private final int nSK;
        private final double[][] ConnAugMatrix;
        private final boolean[] AtomAromatic;
        private final boolean ExplicitHydrogen;
        private double C_026 = 0;
        private double H_051 = 0;
        private double H_052 = 0;


        private GCAtomCentredFragments(IAtomContainer Mol, boolean HasExplicitHydrogen) {

            // Init all variables
            this.CurMol = Mol;
            this.ExplicitHydrogen = HasExplicitHydrogen;
            nSK = Mol.getAtomCount();
            ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
            AtomAromatic = new boolean[nSK];
            for (int i=0; i<nSK; i++)
                AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

            // Check all atoms and try to assign to an ACF
            for (int AtomIdx=0; AtomIdx<nSK; AtomIdx++)
                ProcessAtom(AtomIdx);

        }

        private void ProcessAtom(int AtomIndex) {

            IAtom CurAt = CurMol.getAtom(AtomIndex);

            if (ExplicitHydrogen)
                if (ConnAugMatrix[AtomIndex][AtomIndex] == 1)
                    return;

            if (IsAtomHalogen(CurAt.getAtomicNumber()))
                return;


            int VD=0, nH=0;
            int nSingle=0, nDouble=0, nTriple=0, nArom=0;

            int sX=0, dX=0, tX=0, aX=0, asX=0;
            int sR=0, dR=0, tR=0, aR=0;

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
                        H_051++;

                    } else {

                        // C0sp3 with 1 X atom attached to next C
                        if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==1) && (C_CM==0))
                            H_052++;

                    }
                }
            }

            if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                if ((nH==0) && (sR==4))
                { return; }

                if ((nH==0) && (sR==3) && (sX==1))
                {return; }

                if ((nH==0) && (nDouble==1) && (sR==2) && (dX==0))
                { return; }

                if ((nH==0) && (nDouble==1) && (sX==1) && (sR==1) && (dX==0))
                { return; }

                if ( ((nH==0) && (tR==1) && (sR==1)) ||
                        ((nH==0) && (dR==2)) )
                { return; }

                if ((nH==0) && (tR==1) && (sX==1))
                { return; }

                if ((nH==0) && (VD==3) && (aR>=2) && ((sR + dR ==1)||(aR==3)))
                { return; }

                if ((nH==0) && (VD==3) && (aR==2) && (sX + dX == 1))
                { C_026++; }

            }
        }



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

        private boolean IsAtomHalogen(int AtomicNumber) {
            if ((AtomicNumber==9)||(AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53))
                return true;
            return false;
        }

    }
    private class MoleculePaths {

        // Limit for the number of atoms
        private static final int MAX_PATH_LENGTH = 2000;
        private static final int MAX_PATH_LENGTH_FOR_WALKS = 10;

        private IAtomContainer m;

        // Path indices
        public double[] Path_Counts;
        public double[] Multiple_Path_Counts;
        public double[] TotalPC;
        public double Total_Path_Count;
//        public double IDpi;
//        public double PCR;
//        public double PCD;
//        public double ID_Randic;
//        public double ID_Balaban;
        public double[] Pws;


        // private vars accessible to all methods for iterative DFS
        private double[] Vertex_Distance_Degree;
        private boolean[] Entered;
        private double[] TotPC;
        private double[] TotPCMult;
        private double[][] AdjConnectionMatrix;
        private int nSK;

        private MoleculePaths(InsilicoMolecule Mol) throws GenericFailureException {
            Calculate(Mol);
        }

        private void Calculate(InsilicoMolecule Mol) throws GenericFailureException {

            //// Get structure, matrices and initial settings

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

            Path_Counts = new double[MAX_PATH_LENGTH];
            Multiple_Path_Counts = new double[MAX_PATH_LENGTH];

            double[][] AdjPow = new double[nSK][nSK];

            double[][] AWC = new double[nSK][MAX_PATH_LENGTH];

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
                    int aw = 0;
                    for (int j=0; j< nSK; j++)
                        aw += AdjPow[i][j];
                    AWC[i][Pow_k] = aw;
                }

                Pow_k++;

            }


            //// Paths calculation

            int LenChi = 5;

            int lim = Math.min(nSK, MAX_PATH_LENGTH);
            TotPC = new double[lim];
            TotalPC = new double[lim];
            TotPCMult = new double[lim];
            for (int z=0; z<lim;z++)
                TotPC[z] = 0;

            Vertex_Distance_Degree = new double[nSK];

            for (int i=0; i<nSK; i++) {
                Vertex_Distance_Degree[i] = 0;
                for (int j = 0; j < nSK; j++)
                    Vertex_Distance_Degree[i] += (double) TopoDistMatrix[i][j];
            }

            Entered = new boolean[nSK];
            int jlim = LenChi;
            Pws = new double[jlim];
            for (int i=0; i<jlim; i++)
                Pws[i] = 0;
            if (nSK <= LenChi)
                jlim = nSK - 1;
            double[] TotPCprev = new double[jlim];

            int ii = 0;

            for (int i=0; i<nSK; i++) {

                for (int j = 0; j < nSK; j++)
                    Entered[j] = false;
                Entered[i] = true;

                int PathLength = 0;

//                double VD1 = 0;
//                for (int k=0; k<nSK; k++)
//                    if (k != i)
//                        if (AdjMat[i][k] != 0)
//                            VD1++;

                // Cycle on all atoms connected to i-th
                for (int j = 0; j < nSK; j++) {
                    if (j == i) continue;
                    if (AdjConnectionMatrix[i][j] != 0) {

//                        double VD2 = 0;
//                        for (int k=0; k<nSK; k++)
//                            if (k != j)
//                                if (AdjMat[j][k] != 0)
//                                    VD2++;

                        double Cur_Mult_Bond_Order = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(i), m.getAtom(j)));
//                        double Cur_Mult_Ver_Deg = 1.0 / Math.sqrt(VD1 * VD2);
//                        double Cur_Balaban_Weight = 1.0 / Math.sqrt(Vertex_Distance_Degree[i] * Vertex_Distance_Degree[j]);

//                        NextPathVisit(j, PathLength, Cur_Mult_Bond_Order, Cur_Mult_Ver_Deg, Cur_Balaban_Weight);
                        NextPathVisit(j, PathLength, Cur_Mult_Bond_Order);

                    }
                }

                for (int k=0; k<jlim; k++) {
                    if (AWC[ii][k] != 0)
                        Pws[k] += (TotPC[k] - TotPCprev[k]) / AWC[ii][k];
                    TotPCprev[k] = TotPC[k];
                }

                ii++;
            }

            for (int k=0; k<jlim; k++)
                Pws[k] = Pws[k] / nSK;

            for (int i=0; i<lim; i++) {
                TotPC[i] /= 2.0;
                TotalPC[i] = TotPC[i];
                TotPCMult[i] /= 2.0;
            }

            Total_Path_Count = nSK;
//            IDpi = nSK;
//            ID_Randic = nSK;
//            ID_Balaban = nSK;

            for (int i=0; i<lim; i++) {
                Total_Path_Count = Total_Path_Count + TotPC[i];
//                IDpi = IDpi + TotPCMult[i];
            }

//            PCR = (IDpi / Total_Path_Count);
//            PCD = (IDpi - Total_Path_Count);

            for (int i=0; i<lim; i++) {
                Path_Counts[i] = TotPC[i];
                Multiple_Path_Counts[i] = TotPCMult[i];
            }
        }

//        private void NextPathVisit(int Atom_Idx, int PathLength, double Mult_Bond_Order,
//                                   double Mult_Ver_Deg, double Balaban_Weight) {
        private void NextPathVisit(int Atom_Idx, int PathLength, double Mult_Bond_Order) {

            Entered[Atom_Idx] = true;
            if (PathLength < TotPC.length) {

                TotPC[PathLength] += 1;
                TotPCMult[PathLength] += Mult_Bond_Order;

                for (int Next_Atom_Idx = 0; Next_Atom_Idx < nSK; Next_Atom_Idx++) {
                    if (Next_Atom_Idx == Atom_Idx) continue;
                    if (Entered[Next_Atom_Idx]) continue;
                    if (AdjConnectionMatrix[Atom_Idx][Next_Atom_Idx] != 0) {

//                        double BW1 = Vertex_Distance_Degree[Atom_Idx];
//                        double BW2 = Vertex_Distance_Degree[Next_Atom_Idx];

                        double BO = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(Atom_Idx), m.getAtom(Next_Atom_Idx)));

//                        double VD = 0;
//                        for (int k=0; k<nSK; k++)
//                            if (k != Atom_Idx)
//                                if (AdjConnectionMatrix[Atom_Idx][k] != 0)
//                                    VD++;
//
//                        double VD2= 0;
//                        for (int k=0; k<nSK; k++)
//                            if (k != Next_Atom_Idx)
//                                if (AdjConnectionMatrix[Next_Atom_Idx][k] != 0)
//                                    VD2++;

                        double Cur_Mult_Bond_Order = BO * Mult_Bond_Order;
//                        double Cur_Mult_Ver_Deg = Mult_Ver_Deg * 1.0 / Math.sqrt(VD * VD2);
//                        double Cur_Balaban_Weight = Balaban_Weight * 1.0 / Math.sqrt(BW1 * BW2);

//                        NextPathVisit(Next_Atom_Idx, PathLength + 1, Cur_Mult_Bond_Order, Cur_Mult_Ver_Deg, Cur_Balaban_Weight);
                        NextPathVisit(Next_Atom_Idx, PathLength + 1, Cur_Mult_Bond_Order);

                    }
                }

            }
            Entered[Atom_Idx] = false;

        }

    }
}
