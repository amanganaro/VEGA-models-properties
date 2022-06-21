package insilico.logk.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.basic.WeightsPolarizability;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsQuantumNumber;
import insilico.core.descriptor.blocks.weights.other.WeightsValenceVertexDegree;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
//import insilico.descriptor.blocks.*;

//import insilico.descriptor.blocks.logP.weights.GCWeights;
import insilico.logk.descriptors.weights.GCAtomCentredFragments;
import insilico.logk.descriptors.weights.GCWeights;
import insilico.logk.descriptors.weights.MoleculePaths;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

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

    private double MISSING_VALUE = -999;

    public double ALogP = MISSING_VALUE;
    public double P_VSA_i_2 = MISSING_VALUE;
    public double MLogP = MISSING_VALUE;
    public double P_VSA_p_3 = MISSING_VALUE;
    public double Eta_betaP = MISSING_VALUE;
    public double CATS2D_00_LL = MISSING_VALUE;
    public double C = MISSING_VALUE;
    public double PCD = MISSING_VALUE;
    public double Ui = MISSING_VALUE;
    public double N = MISSING_VALUE;
    public double totalcharge = MISSING_VALUE;
    public double CATS2D_00_PP = MISSING_VALUE;
    public double C_024 = MISSING_VALUE;
    public double nCsp2 = MISSING_VALUE;

    public double MLOGP2 = MISSING_VALUE;
    public double GATS1i = MISSING_VALUE;
    public double SpMax2_Bh_P = MISSING_VALUE;
    public double nBm = MISSING_VALUE;
    public double MATS5e = MISSING_VALUE;
    public double AMW = MISSING_VALUE;
    public double F01_C_N = MISSING_VALUE;
    public double T_O_O = MISSING_VALUE;
    public double J_D_DT = MISSING_VALUE;
    public double SpMax_AEA_dm = MISSING_VALUE;

    private static final int MAX_PATH_LENGTH = 2000;
    private static final int MAX_PATH_LENGTH_FOR_WALKS = 10;

    private IAtomContainer m;

    // Walk indices
    public double[] Walk_Counts;
    public double[] Self_Returning_Walk_Counts;
    public double Total_Walk_Count;

    // Path indices
    public double[] Path_Counts;
    public double[] Multiple_Path_Counts;
    public double[] TotalPC;
    public double Total_Path_Count;
    public double IDpi;
    public double PCR;
    public double ID_Randic;
    public double ID_Balaban;
    public double[] Pws;


    // private vars accessible to all methods for iterative DFS
    private double[] Vertex_Distance_Degree;
    private boolean[] Entered;
    private double[] TotPC;
    private double[] TotPCMult;
    private double[] IDRandic;
    private double[] IDBalaban;
    private double[][] AdjConnectionMatrix;

    public double ExperimentalValue;

    public double[] getDescriptors(){
        return new double[] {
                ALogP, P_VSA_i_2, MLogP, P_VSA_p_3, Eta_betaP, CATS2D_00_LL, C,
                PCD, Ui, N, totalcharge, CATS2D_00_PP, C_024, nCsp2, MLOGP2, GATS1i, SpMax2_Bh_P,
                nBm, MATS5e, AMW, F01_C_N, T_O_O, J_D_DT, SpMax_AEA_dm};
    }


    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateLogP(mol);
        CalculateP_VSA(mol);
        CalculateEtaBeta(mol);
        CalculateCATS2D(mol);
        CalculateConstitutional(mol);
        CalculateWAP(mol);
        CalculateACF(mol);
        CalculateAutoCorrelation(mol);
        CalculateSP(mol);
        CalculateSPMax(mol);
        CalculateAtomPairs(mol);
        CalculateJDDT(mol);
        CalculateNcs(mol);
    }

    private void CalculateNcs(InsilicoMolecule mol) {
        String ncsSmart = "[$([C;D2]([#6])[#6]),$([C;D3]([#6])([#6])[*;!#6]),$([C;D4]([#6])([#6])([*;!#6])[*;!#6])]";
        Pattern query = SmartsPattern.create(ncsSmart).setPrepare(false);
        boolean status;
        try {
            status = query.matches(mol.GetStructure());
            if (status) {
                nCsp2 = query.matchAll(mol.GetStructure()).countUnique();
            } else {
                nCsp2 = 0;
            }
        } catch (InvalidMoleculeException ex){
            log.warn(ex.getMessage());
        }


    }

    private void CalculateJDDT(InsilicoMolecule mol) {


        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            J_D_DT = MISSING_VALUE;
            return;
        }
        // Adj matrix available for calculations
        int[][] AdjMatrix = null;
        try {
            AdjMatrix = mol.GetMatrixAdjacency();
        } catch (GenericFailureException ex){
            log.warn(ex.getMessage());
        }

        int nSK = curMol.getAtomCount();
        int nBo = curMol.getBondCount();
        double[][] Mat = null;
        try {
            Mat = mol.GetMatrixDistanceDetour();
        } catch (GenericFailureException ex){
            log.warn(ex.getMessage());
        }

        Matrix DataMatrix = new Matrix(Mat);
        double[] eigenvalues = null;
        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            log.warn("Unable to calculate eigenvalues: " + e.getMessage());
        }

        // Eigenvalue based descriptors
        J_D_DT = 0;
        for (double val : eigenvalues) {
            J_D_DT += Math.abs(val);
        }
    }


    private void CalculateAtomPairs(InsilicoMolecule mol) {

        int MAX_TOPO_DISTANCE = 10;
        String[][] atomCoupleCN = {{"C", "N"}};
        String[][] atomCoupleOO = {{"O", "O"}};

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            F01_C_N = MISSING_VALUE;
            T_O_O = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn("Invalid structure, unable to calculate: " + mol.GetSMILES());
            F01_C_N = MISSING_VALUE;
            T_O_O = MISSING_VALUE;
            return;
        }

        for (int d = 0; d< atomCoupleCN.length; d++) {

            int descT = 0;
            int[] descB = new int[MAX_TOPO_DISTANCE];
            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCoupleCN[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCoupleCN[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2) // DA VEDERE PERCHE MAGGIORE DI 2
                                descT += TopoMat[i][j];

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j]-1] = 1;
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (atomCoupleCN[d][0].compareTo(atomCoupleCN[d][1]) == 0) {
                descT /= 2;
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            F01_C_N = descF[0];
        }

        for (int d = 0; d< atomCoupleOO.length; d++) {

            int descT = 0;
            int[] descB = new int[MAX_TOPO_DISTANCE];
            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCoupleOO[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCoupleOO[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2) // DA VEDERE PERCHE MAGGIORE DI 2
                                descT += TopoMat[i][j];

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j]-1] = 1;
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (atomCoupleOO[d][0].compareTo(atomCoupleOO[d][1]) == 0) {
                descT /= 2;
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            T_O_O = descT;
        }
    }

    private void CalculateSPMax(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            SpMax_AEA_dm = MISSING_VALUE;
            return;
        }

        int nSK = curMol.getAtomCount();
        int nBO = curMol.getBondCount();

        // Only for mol with nSK>1
        if (nSK < 2) {
            log.warn("Less then 2 atoms");
            SpMax_AEA_dm = MISSING_VALUE;
            return;
        }

        // Gets basic matrix
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SpMax_AEA_dm = MISSING_VALUE;
            return;
        }

        double[][] curDataMatrix = new double[nBO][nBO];
        // plain EA
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                curDataMatrix[i][j] = EdgeAdjMat[i][j][0];

        double[] w = new double[nBO];
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

        Matrix DataMatrix = new Matrix(curDataMatrix);
        double[] eigenvalues;
        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            log.warn("Unable to calculate eigenvalues: " + e.getMessage());
            return;
        }

        // check on precision
        for (int i=0; i<eigenvalues.length; i++) {
            if (Math.abs(eigenvalues[i]) < 0.000001)
                eigenvalues[i] = 0;
        }

        // Eigenvalue based descriptors
        double EigMax = eigenvalues[0];
        double EigMin = eigenvalues[0];

        for (double val : eigenvalues) {
//            EigAve += val;
            if (val > EigMax)
                EigMax = val;
            if (val< EigMin)
                EigMin = val;
        }
        double NormEigMax = EigMax / (double) nBO;


        SpMax_AEA_dm = NormEigMax;

    }

    private void CalculateSP(InsilicoMolecule mol) {

        // Burden matrix is calculated on H-filled molecules
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            SpMax2_Bh_P = MISSING_VALUE;
            return;
        }

        // Gets matrix
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            SpMax2_Bh_P = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();
        WeightsPolarizability weightsPolarizability = new WeightsPolarizability();
        double[] w = weightsPolarizability.getScaledWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            if(i != 1){
                continue;
            } else {
                double valH;
                if (i>(eigenvalues.length-1)) {
                    valH = 0;
                } else {
                    if (eigenvalues[eigenvalues.length-1-i] > 0)
                        valH = eigenvalues[eigenvalues.length-1-i];
                    else
                        valH = 0;
                }

                SpMax2_Bh_P = valH;
            }
        }
    }

    private void CalculateAutoCorrelation(InsilicoMolecule mol) {
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            GATS1i = MISSING_VALUE;
            MATS5e = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            GATS1i = MISSING_VALUE;
            MATS5e = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        WeightsIonizationPotential weightsIState = new WeightsIonizationPotential();
        double[] w_i = weightsIState.getScaledWeights(m);


        WeightsElectronegativity weightsElectronegativity = new WeightsElectronegativity();
        double[] w_e = weightsElectronegativity.getScaledWeights(m);

        double wA_i = 0;
        for (int i=0; i<nSK; i++)
            wA_i += w_i[i];
        wA_i = wA_i / ((double) nSK);

        double wA_e = 0;
        for (int i=0; i<nSK; i++)
            wA_e += w_e[i];
        wA_e = wA_e / ((double) nSK);


        // i

        for (int lag=1; lag<=8; lag++) {

            if(lag != 1)
                continue;

            double AC=0, ACS=0, MoranAC=0, GearyAC=0;
            double denom = 0, delta = 0;

            for (int i=0; i<nSK; i++) {

                denom += Math.pow((w_i[i] - wA_i), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        MoranAC += (w_i[i] - wA_i) * (w_i[j] - wA_i);
                        GearyAC += Math.pow((w_i[i] - w_i[j]), 2);
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


            GATS1i = GearyAC;

        }

        for (int lag=1; lag<=8; lag++) {

            if(lag != 5)
                continue;

            double MoranAC=0;
            double denom = 0;
            double delta = 0;

            for (int i=0; i<nSK; i++) {

                denom += Math.pow((w_e[i] - wA_e), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
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
        }

    }

    private void CalculateACF(InsilicoMolecule mol) {
        AtomCenteredFragments block = new AtomCenteredFragments();
        try {
            block.Calculate(mol);
            C_024 = block.GetByName("C-024").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateWAP(InsilicoMolecule mol) {
        try {
            MoleculePaths paths = new MoleculePaths(mol);

            double TPC = paths.Total_Path_Count;
            double piID = paths.IDpi;

            TPC = Math.log(1 + TPC);
            piID = Math.log(1 + piID);

            PCD = piID - TPC;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateConstitutional(InsilicoMolecule mol) {

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nTotH=0;
            int nC=0, nN=0;
            int nHet=0;
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
                if (CurAt.getSymbol().equalsIgnoreCase("N"))
                    nN++;


            }

            N = (nN/(double)(nSK + nTotH))*100;
            C = (nC/(double)(nSK + nTotH))*100;

            //// Counts on bonds

            double scbo=0, nMulBonds = 0;

            for (int i=0; i<nBO; i++) {

                IBond CurBo = curMol.getBond(i);

                if (CurBo.getFlag(CDKConstants.ISAROMATIC)) {
                    nMulBonds++;

                    scbo += 1.5;
                } else {
                    if (CurBo.getOrder() == IBond.Order.SINGLE) {
                        scbo++;
                    } else {
                        nMulBonds++;

                        if (CurBo.getOrder() == IBond.Order.DOUBLE) {
                            scbo += 2;
                        }
                        if (CurBo.getOrder() == IBond.Order.TRIPLE) {
                            scbo += 3;
                        }
                    }
                }

            }

            nBm = nMulBonds;

            WeightsMass weightsMass = new WeightsMass();
            double[] weights = weightsMass.getScaledWeights(curMol);
            double weightH = weightsMass.getScaledWeight("H");

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

            double ave = MISSING_VALUE;
            if (sum != MISSING_VALUE)
                ave = sum/(nSK + nTotH);

            AMW = ave;

            Ui = Math.log(1 + scbo - nBO) / Math.log(2);
            totalcharge = 0;
            for (IAtom at : mol.GetStructure().atoms())
                totalcharge += at.getFormalCharge();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            N = MISSING_VALUE;
            C = MISSING_VALUE;
            nBm = MISSING_VALUE;
            AMW = MISSING_VALUE;
            Ui = MISSING_VALUE;
            totalcharge = MISSING_VALUE;
        }
    }

    private void CalculateCATS2D(InsilicoMolecule mol) {

        String[] TYPE_P = { "P", ""};
        String[] TYPE_L = { "L", ""};

        String[][][] AtomCouples = {
                {TYPE_P, TYPE_P},
                {TYPE_L, TYPE_L}
        };


        IAtomContainer curMol = null;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_00_PP  = MISSING_VALUE;
            CATS2D_00_LL = MISSING_VALUE;
            log.warn("Invalid structure: " + mol.GetSMILES());
            e.printStackTrace();
        }

        int nSK = curMol.getAtomCount();


        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_00_PP  = MISSING_VALUE;
            CATS2D_00_LL = MISSING_VALUE;
            return;
        }

        double[][] ConnAugMatrix = null;
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_00_PP  = MISSING_VALUE;
            CATS2D_00_LL = MISSING_VALUE;
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(curMol, ConnAugMatrix);
        for (String[][] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                for (int j = i; j < nSK; j++) {
                    if (isInCatsType(CatsTypes[i], CatsTypes[j], atomCouple[0][0], atomCouple[1][0]))
                        if (TopoMat[i][j] < 10)
                            desc[TopoMat[i][j]]++;

                }
            }


            for (int i = 0; i < desc.length; i++) {
                if(atomCouple[0][0].equals("L") && i == 0)
                    CATS2D_00_LL = desc[i];
                if(atomCouple[0][0].equals("P") && i == 0)
                    CATS2D_00_PP = desc[i];
            }
        }
    }

    private void CalculateEtaBeta(InsilicoMolecule mol) {

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure, unable to calculate: " + mol.GetSMILES());
            Eta_betaP = MISSING_VALUE;
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
            Eta_betaP = MISSING_VALUE;
            return;
        }

        RingSet MolRings;
        try {
            MolRings = mol.GetSSSR();
        } catch (Exception e) {
            log.warn(e.getMessage());
            Eta_betaP = MISSING_VALUE;
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
                Eta_betaP = MISSING_VALUE;
                return;
            }
            isAromatic[i] = curMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);
        }



        //// Parameters for ETA indices

        double[] alpha = new double[nSK];
        double[] beta_s = new double[nSK];
        double[] beta_ns = new double[nSK];
        double[] gamma = new double[nSK];
        double[] epsilon = new double[nSK];

        for (int i=0; i<nSK; i++) {
            alpha[i] = 0;
            beta_s[i] = 0;
            beta_ns[i] = 0;
            gamma[i] = 0;
            epsilon[i] = 0;
        }

        for (int i=0; i<nSK; i++) {
            alpha[i] = (double)(AtomicNumber[i] - ElectronNumber[i]) / ElectronNumber[i] * (1.0 / (QuantumNumber[i] - 1));
            epsilon[i] = -1.0 * alpha[i] + 0.3 * ElectronNumber[i];
        }

        for (int i=0; i<nSK; i++) {
            double correctionFactor = 0;

            // beta_ns for aromatic structures
            if (isAromatic[i]) {
                // N in pyrrole
                boolean exception = curMol.getAtom(i).getSymbol().equalsIgnoreCase("N") && (VertexDegree[i] == 3);
                // P in pyrrole
                exception = exception || (curMol.getAtom(i).getSymbol().equalsIgnoreCase("P") && (VertexDegree[i] == 3));
                // S, O, Se
                exception = exception || (curMol.getAtom(i).getSymbol().equalsIgnoreCase("O") ||
                        curMol.getAtom(i).getSymbol().equalsIgnoreCase("S") || curMol.getAtom(i).getSymbol().equalsIgnoreCase("Se") );

                if (!exception)
                    beta_ns[i] += 2;
            }

            // beta params
            // cycle on all connected atoms
            for (int at_idx=0; at_idx<nSK; at_idx++) {
                if (at_idx == i) continue;
                if (AdjConnMatrix[i][at_idx] == 0) continue;

                double DeltaEpsilon = Math.abs(epsilon[i] - epsilon[at_idx]);

                // beta_s - consider sigma bonds (all bonds)
                if (DeltaEpsilon <= 0.3)
                    beta_s[i] += 0.5;
                else
                    beta_s[i] += 0.75;

                // beta_ns - consider dbl and triple pi bonds

                // correction factor (O, S with tot bnd order = 2, N, P with tot bond order = 3
                if ( (!isAromatic[i]) && (isAromatic[at_idx]) ) {
                    int TotBondOrder = 0;
                    for (int k=0; k<nSK; k++) {
                        if (k == i) continue;
                        TotBondOrder += AdjConnMatrix[i][k];
                    }
                    // tot bond order includes H??
                    TotBondOrder += Manipulator.CountImplicitHydrogens(curMol.getAtom(i));

                    if ( (curMol.getAtom(i).getSymbol().equalsIgnoreCase("O") ) || (curMol.getAtom(i).getSymbol().equalsIgnoreCase("S") ) )
                        if (TotBondOrder == 2)
                            correctionFactor = 0.5;
                    if ( (curMol.getAtom(i).getSymbol().equalsIgnoreCase("N") ) || (curMol.getAtom(i).getSymbol().equalsIgnoreCase("P") ) )
                        if (TotBondOrder == 3)
                            correctionFactor = 0.5;

                }

                // bonds contribution - triple bnd count as 1 sigma and 2 pi
                if (AdjConnMatrix[i][at_idx] == 3) {
                    if (DeltaEpsilon <= 0.3)
                        beta_ns[i] += (1 + 1); // similar electronegativity
                    else
                        beta_ns[i] += (1.5 + 1.5); // different electronegativity
                }

                // bonds contribution - dbl bnd count as 1 sigma and 1 pi
                if (AdjConnMatrix[i][at_idx] == 2) {

                    // check conjugation
                    boolean SameAromRing = false;

                    for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
                        IRing curRing = (IRing) MolRings.getAtomContainer(r);
                        if ( (curRing.contains(curMol.getAtom(i))) && (curRing.contains(curMol.getAtom(at_idx))) ) {
                            if ( (curMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) && (curMol.getAtom(at_idx).getFlag(CDKConstants.ISAROMATIC)) ) {
                                SameAromRing = true;
                                break;
                            }
                        }
                    }

                    if ( (!SameAromRing) && (isConjugated(curMol, i, AdjConnMatrix)) && (isConjugated(curMol, at_idx, AdjConnMatrix)) ) {

                        beta_ns[i] += (1.5);

                    } else {

                        if (DeltaEpsilon <= 0.3)
                            beta_ns[i] += (1); // similar electronegativity
                        else
                            beta_ns[i] += (1.5); // different electronegativity
                    }
                }
            }

            // Correction for NO2 for beta_ns
            if (AdjConnMatrix[i][i] == 7) {
                if (curMol.getAtom(i).getFormalCharge() == +1) {
                    int Odbl = 0, Osngminus = 0, vd = 0;
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (AdjConnMatrix[i][j] != 0) {
                            vd++;
                            if (AdjConnMatrix[j][j] == 8) {
                                if (AdjConnMatrix[i][j] == 2)
                                    Odbl++;
                                if (AdjConnMatrix[i][j] == 1)
                                    if (curMol.getAtom(j).getFormalCharge() == -1)
                                        Osngminus++;
                            }
                        }
                    }

                    if ((vd == 3) && (Osngminus == 1) && (Odbl == 1)) {
                        // NO2 in O=[N+]-[O-] form
                        // for compliance with Dragon 7, a correction is used as this group would be
                        // considered as a O=N=O
                        beta_ns[i] += 1.5 + 1.5;
                    }
                }
            }


            beta_ns[i] += correctionFactor;
//            gamma[i] = alpha[i] / (beta_s[i] + beta_ns[i]);
//            beta_s[i] /= 2.0;
            beta_ns[i] = (beta_ns[i] - correctionFactor) / 2.0 + correctionFactor;

        }

        //// ETA indices

        double eta_betans = 0;

        for (int i=0; i<nSK; i++) {
            eta_betans += beta_ns[i];
        }

        Eta_betaP = eta_betans;
    }




    private void CalculateP_VSA(InsilicoMolecule mol) {
//        DescriptorBlock block = new P_VSA();

        String[][] PPP_TYPES = {
                {"D",""},
                {"A",""},
                {"P",""},
                {"N",""},
                {"L",""},
                {"Cyc",""},

        };

        IAtomContainer curMolNoH;
        double[][] ConnAugMatrixNoH = null;
        try {
            curMolNoH = mol.GetStructure();
            ConnAugMatrixNoH = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
            P_VSA_p_3= MISSING_VALUE;
            P_VSA_i_2 = MISSING_VALUE;
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
            log.warn("Invalid structure: " + mol.GetSMILES());
            P_VSA_p_3= MISSING_VALUE;
            P_VSA_i_2 = MISSING_VALUE;
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
            int at2 = m.getAtomNumber(bo.getAtom(1));
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

        WeightsIonizationPotential weightsIonizationPotential = new WeightsIonizationPotential();
        WeightsPolarizability weightsPolarizability = new WeightsPolarizability();
        double[] w_i = weightsIonizationPotential.getScaledWeights(m);
        double[] w_p = weightsPolarizability.getScaledWeights(m);

        // i // p
        int bins = 4;
        for (int b=1; b<=bins; b++) {
            double PVSA = 0;
            for (int i = 0; i < nSK; i++) {
                if (w_i[i] == Descriptor.MISSING_VALUE) continue;
                if (b == CalculateBin("i", w_i[i]))
                    PVSA += VSA[i];
            }
            if(b == 2)
                P_VSA_i_2 = PVSA;
        }
        for (int b=1; b<=bins; b++) {
            double PVSA = 0;
            for (int i = 0; i < nSK; i++) {
                if (w_p[i] == Descriptor.MISSING_VALUE) continue;
                if (b == CalculateBin("p", w_p[i]))
                    PVSA += VSA[i];
            }
            if(b == 3)
                P_VSA_p_3 = PVSA;
        }
    }

    private void CalculateLogP(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ALogP = MISSING_VALUE;
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF = new AtomCenteredFragments();
            CurACF.Calculate(mol);


        double LogP = 0;
        double MR = 0;
        int[] Frags;
        GCAtomCentredFragments GCfrags;

        try {
            GCfrags = new GCAtomCentredFragments(mol.GetStructure(), false);
            Frags = GCfrags.GetACF();
        } catch (InvalidMoleculeException e) {
            ALogP = MISSING_VALUE;
            return;
        }

        // Check if some fragments are missing values
        for (int f : Frags)
            if (f == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += GCWeights.GetHydrophobiticty(Frags[i]);
            MR += GCWeights.GetMolarRefractivity(Frags[i]);
        }

        for (int k : GCfrags.getNotMappedFragCount().keySet()) {
            LogP += GCWeights.GetHydrophobiticty(k) * GCfrags.getNotMappedFragCount().get(k);
            MR += GCWeights.GetMolarRefractivity(k) * GCfrags.getNotMappedFragCount().get(k);

        }

        ALogP = LogP;

        MLOGP2 = Math.pow(MLogP, 2);
    }

    private int CalculateBin(String curWeight, double value) {



        if (curWeight.equalsIgnoreCase("p")) {
            // Polarizability
            if (value < 0.4) return 1;
            if (value < 1.0) return 2;
            if (value < 2.0) return 3;
            return 4;
        }

        if (curWeight.equalsIgnoreCase("i")) {
            // Ionization potential
            if (value < 1) return 1;
            if (value < 1.15) return 2;
            if (value < 1.25) return 3;
            return 4;
        }


        return 0;
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

    private void NextPathVisit(int Atom_Idx, int PathLength, double Mult_Bond_Order,
                               double Mult_Ver_Deg, double Balaban_Weight, IAtomContainer m) {

        Entered[Atom_Idx] = true;
        int nSK = m.getAtomCount();
        if (PathLength < TotPC.length) {

            TotPC[PathLength] += 1;
            TotPCMult[PathLength] += Mult_Bond_Order;
            IDRandic[PathLength] += Mult_Ver_Deg;
            IDBalaban[PathLength] += Balaban_Weight;

            for (int Next_Atom_Idx = 0; Next_Atom_Idx < nSK; Next_Atom_Idx++) {
                if (Next_Atom_Idx == Atom_Idx) continue;
                if (Entered[Next_Atom_Idx]) continue;
                if (AdjConnectionMatrix[Atom_Idx][Next_Atom_Idx] != 0) {

                    double BW1 = Vertex_Distance_Degree[Atom_Idx];
                    double BW2 = Vertex_Distance_Degree[Next_Atom_Idx];

                    double BO = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(Atom_Idx), m.getAtom(Next_Atom_Idx)));

                    double VD = 0;
                    for (int k=0; k<nSK; k++)
                        if (k != Atom_Idx)
                            if (AdjConnectionMatrix[Atom_Idx][k] != 0)
                                VD++;

                    double VD2= 0;
                    for (int k=0; k<nSK; k++)
                        if (k != Next_Atom_Idx)
                            if (AdjConnectionMatrix[Next_Atom_Idx][k] != 0)
                                VD2++;



                    double Cur_Mult_Bond_Order = BO * Mult_Bond_Order;
                    double Cur_Mult_Ver_Deg = Mult_Ver_Deg * 1.0 / Math.sqrt(VD * VD2);
                    double Cur_Balaban_Weight = Balaban_Weight * 1.0 / Math.sqrt(BW1 * BW2);

                    NextPathVisit(Next_Atom_Idx, PathLength + 1, Cur_Mult_Bond_Order, Cur_Mult_Ver_Deg, Cur_Balaban_Weight , m);

                }
            }

        }
        Entered[Atom_Idx] = false;

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
                log.warn("unable to count H");
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



}
