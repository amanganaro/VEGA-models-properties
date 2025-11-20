package insilico.dio1.model.modelGBTpmml;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.basic.WeightsPolarizability;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.*;

public class DioDescriptors {

    public double[] Calculate(String SMILES) throws GenericFailureException {
        InsilicoMolecule mol = SmilesMolecule.Convert(SMILES);
        return Calculate(mol);
    }

    public double[] Calculate(InsilicoMolecule mol) throws GenericFailureException {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            throw new GenericFailureException(e);
        }
        int nSK = m.getAtomCount();
        int nBO = m.getBondCount();

        int[][] TopoMatrix;
        double[][] BurdenMatrix;
        double[][] DistDetMatrix;
        int[][] LaplaceMatrix;
        try {
            TopoMatrix =  mol.GetMatrixTopologicalDistance();
            BurdenMatrix = mol.GetMatrixBurden();
            DistDetMatrix = mol.GetMatrixDistanceDetour();
            LaplaceMatrix = mol.GetMatrixLaplace();
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate matrices");
        }

        // I and E states
        WeightsIState w_istate = new WeightsIState();
        double[] w_is = w_istate.getWeights(m, false);
        for (double val : w_is)
            if (val == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("unable to calculate all E-states");
        double[] w_es = new double[nSK];
        for (int at = 0; at<nSK; at++) {
            double sumDeltaI = 0;
            for (int j = 0; j < nSK; j++)
                if (at != j)
                    sumDeltaI += (w_is[at] - w_is[j]) / Math.pow((double) TopoMatrix[at][j] + 1.0, 2.0);
            w_es[at] = w_is[at] + sumDeltaI;
        }



        //// TI2_L

        double[][] mat = new double[LaplaceMatrix.length][LaplaceMatrix[0].length];
        for (int i=0; i<LaplaceMatrix.length; i++)
            for (int j=0; j<LaplaceMatrix[0].length; j++)
                mat[i][j] = LaplaceMatrix[i][j];
        double D02_TI2_L = CalculateEigDescriptors(mat, nSK).get("TI2");



        //// D/Dt

        double D03_SpPosA_D_Dt = CalculateEigDescriptors(DistDetMatrix, nSK).get("SpPosA");


        //// AVS + HyWI on B(e)

        double[] w = (new WeightsElectronegativity()).getScaledWeights(m) ;
        for (int i=0; i<nSK; i++)
            BurdenMatrix[i][i] = w[i];

        double AVS = 0;
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++)
                AVS += BurdenMatrix[i][j];
        double D04_HAVS_B_e = AVS / (double) nSK;

        double hywi = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i; j<nSK; j++)
                hywi += ( Math.pow(BurdenMatrix[i][j], 2) + BurdenMatrix[i][j]);
        double D05_HyWi_B_e = Math.log(1 + 0.5 * hywi);


        // EE on B(s)

        for (int i=0; i<nSK; i++)
            BurdenMatrix[i][i] = w_is[i];
        double D06_EE_B_s = CalculateEigDescriptors(BurdenMatrix, nSK).get("EE");


        //// SsCH3

        double D17_SsCH3 = 0;
        for (int at=0; at<m.getAtomCount(); at++) {
            IAtom curAt = m.getAtom(at);
            if (curAt.getAtomicNumber() != 6)
                continue;

            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) { }

            int Charge;
            try {
                Charge = curAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }

            // Count bonds
            int nBnd=0, nSng = 0, nAr=0;
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
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nBnd++;
                }
            }

            // aaCH group
            if ( (nBnd==1) && (nSng==1) && (nAr == 0) && (nH == 3) && (Charge == 0) )
                D17_SsCH3 += w_es[at];

        }



        //////// descriptors with H-filled structure

        IAtomContainer curMolH = Manipulator.AddHydrogens(m);

        int[][] TopoDistMatH;
        double[][] ConnMatH;
        double[][] ConnAugMatH;
        try {
            ConnMatH= ConnectionAugMatrix.getMatrix(curMolH);
            ConnAugMatH = ConnectionAugMatrix.getMatrix(curMolH);
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate matrix - " + e);
        }
        int nSKH = curMolH.getAtomCount();
        int nBOH = curMolH.getBondCount();
        int[] VertexDegH = WeightsVertexDegree.getWeights(curMolH, false);


        //// SIC3

        int CurLag = 3;
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

        double ic=0;
        for (int i=0; i<Gn.size(); i++)
            ic += ((double)Gn.get(i)/nSKH) * (Math.log((double)Gn.get(i)/nSKH));
        ic = (-1.00 / Math.log(2)) * ic;
        double diff = Math.log(nSKH) / Math.log(2);
        double sic = ic / diff;

        double D01_SIC3 = sic;


        //// Autocorrelations

        w = (new WeightsMass()).getScaledWeights(curMolH);
        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing mass weight");
        HashMap<String, Double> auto_m = Autocorrelation(w, TopoDistMatH, nSKH);

        w = (new WeightsElectronegativity()).getScaledWeights(curMolH);
        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing electronegativity weight");
        HashMap<String, Double> auto_e = Autocorrelation(w, TopoDistMatH, nSKH);

        w = (new WeightsPolarizability()).getScaledWeights(curMolH);
        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing polarizability weight");
        HashMap<String, Double> auto_p = Autocorrelation(w, TopoDistMatH, nSKH);

        w = (new WeightsIonizationPotential()).getScaledWeights(curMolH);
        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing ionization potential weight");
        HashMap<String, Double> auto_i = Autocorrelation(w, TopoDistMatH, nSKH);

        w = (new WeightsIState()).getWeights(curMolH, true);
        for (int i=0; i<nSKH; i++) {
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing I-State weight");
            if (curMolH.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }
        HashMap<String, Double> auto_s = Autocorrelation(w, TopoDistMatH, nSKH);

        double D07_ATSC2m = auto_m.get("ATSC2");
        double D08_MATS7m = auto_m.get("MATS7");
        double D09_MATS2e = auto_e.get("MATS2");
        double D10_MATS2i = auto_i.get("MATS2");
        double D11_MATS1s = auto_s.get("MATS1");
        double D12_GATS3e = auto_e.get("GATS3");
        double D13_GATS8e = auto_e.get("GATS8");
        double D14_GATS5p = auto_p.get("GATS5");
        double D15_GATS6i = auto_i.get("GATS6");


        //// P_VSA

        double[] VSA = new double[nSKH];
        double[] PartialVSA = new double[nSKH];
        double[] IdealBondLen = new double[nBOH];
        double[] VdWRadius = new double[nSKH];

        for (int i=0; i<nSKH; i++) {
            IAtom at = curMolH.getAtom(i);
            VSA[i] = 0;
            PartialVSA[i] = 0;
            VdWRadius[i] = GetVdWRadius(curMolH, at);
            if (VdWRadius[i] == Descriptor.MISSING_VALUE) {
                VSA[i] = Descriptor.MISSING_VALUE;
                PartialVSA[i] = Descriptor.MISSING_VALUE;
            }
        }

        for (int i=0; i<nBOH; i++) {
            IBond bo = curMolH.getBond(i);
            double refR = this.GetRefBondLength(bo.getAtom(0), bo.getAtom(1));

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

        for (int b=0; b<nBOH; b++) {
            IBond bo = curMolH.getBond(b);

            int at1 = curMolH.indexOf(bo.getAtom(0));
            int at2 = curMolH.indexOf(bo.getAtom(1));
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

        for (int i=0; i<nSKH; i++) {
            if (VSA[i] != Descriptor.MISSING_VALUE)
                VSA[i] = 4 * Math.PI * Math.pow(VdWRadius[i], 2) - (Math.PI * VdWRadius[i] * PartialVSA[i]);
            if (VSA[i] < 0)
                VSA[i] = 0;
        }

        // in the new mol the H are all added after the original atoms, so the first nSK atoms
        // and calculated VSA are the same both for the original mol and for the H-filled one

        double D16_P_VSA_ppp_L = 0;
        for (int i = 0; i < nSK; i++) {
            if (isLtype(i, curMolH, ConnAugMatH))
                D16_P_VSA_ppp_L += VSA[i];
        }



        double[] res = new double[17];
        res[0] = D01_SIC3;
        res[1] = D02_TI2_L;
        res[2] = D03_SpPosA_D_Dt;
        res[3] = D04_HAVS_B_e;
        res[4] = D05_HyWi_B_e;
        res[5] = D06_EE_B_s;
        res[6] = D07_ATSC2m;
        res[7] = D08_MATS7m;
        res[8] = D09_MATS2e;
        res[9] = D10_MATS2i;
        res[10] = D11_MATS1s;
        res[11] = D12_GATS3e;
        res[12] = D13_GATS8e;
        res[13] = D14_GATS5p;
        res[14] = D15_GATS6i;
        res[15] = D16_P_VSA_ppp_L;
        res[16] = D17_SsCH3;

        return res;
    }



    public final static String[] DESCRIPTOR_NAMES = {
        "SIC3",
        "TI2_L",
        "SpPosA_D/Dt",
        "AVS_B(e)",
        "HyWi_B(e)",
        "EE_B(s)",
        "ATSC2m",
        "MATS7m",
        "MATS2e",
        "MATS2i",
        "MATS1s",
        "GATS3e",
        "GATS8e",
        "GATS5p",
        "GATS6i",
        "P_VSA_ppp_L",
        "SsCH3"
    };

    // values calculated from the saved TRAINING SET only
    private final double[] DESC_AVE =
            {0.785820, 3.187552, 0.765788, 3.375498, 3.595550, 7.539290, 13.572735, -0.011555, 0.111526, -0.019826, -0.119320, 1.013353, 0.744149, 0.926210, 0.790203, 70.611304, 3.116276};
    private final double[] DESC_STD =
            {0.111934, 1.772228, 0.143177, 0.237223, 0.502346, 1.194815, 9.091431, 0.445384, 0.196965, 0.198687, 0.114697, 0.305882, 0.774568, 0.357813, 0.411515, 51.039355, 3.280983};

    public double[] Scale(double[] descriptors) {
        double[] scaled = new double[descriptors.length];
        for (int i=0; i<descriptors.length; i++)
            scaled[i] = ( descriptors[i] - DESC_AVE[i] ) / DESC_STD[i];
        return scaled;
    }



    private HashMap<String, Double> CalculateEigDescriptors(double[][] mat, int nSK) throws GenericFailureException {
        HashMap<String, Double> res = new HashMap<>();

        // SpPosA
        // EE

        Matrix DataMatrix = new Matrix(mat);
        double[] eigenvalues;
        Matrix eigenvectors;
        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            throw new GenericFailureException("unable to calculate eigenvalues - " + e);
        }

        // check on precision
        for (int i=0; i<eigenvalues.length; i++) {
            if (Math.abs(eigenvalues[i]) < 0.000001)
                eigenvalues[i] = 0;
        }

        double LastEigen = 0;
        int matrix_dim = eigenvalues.length;
        for (int i=eigenvalues.length-1; i>=0; i--)
            if (eigenvalues[i] > 0)
                LastEigen = eigenvalues[i];
        double TI2;
        if (LastEigen > 0)
            TI2 = 4.0 / (matrix_dim * LastEigen);
        else
            TI2 = 0;
        res.put("TI2", TI2);

        double SpPos = 0, EigExpSum = 0;
        for (double val : eigenvalues) {
            if (val > 0)
                SpPos += val;
            EigExpSum += Math.exp(val);
        }
        double SpPosA = SpPos / (double) nSK;
        res.put("SpPosA", SpPosA);
        double EstradaLike = Math.log(1 + EigExpSum);
        res.put("EE", EstradaLike);

        return res;
    }


    private HashMap<String, Double> Autocorrelation(double[] w, int[][] TopoMatrix, int nSKH) {
        HashMap<String, Double> res = new HashMap<>();

        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSKH; i++)
            wA += w[i];
        wA = wA / ((double) nSKH);

        for (int lag=1; lag<=8; lag++) {

            double ACS=0, MoranAC=0, GearyAC=0;
            double denom = 0, delta = 0;

            for (int i=0; i<nSKH; i++) {

                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSKH; j++)
                    if (TopoMatrix[i][j] == lag) {
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
                    MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSKH)) * denom);
                    GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSKH - 1))) * denom);
                }
            }

            ACS /= 2.0;

            // Sets descriptors
            res.put("ATSC" + lag, ACS);
            res.put("MATS" + lag, MoranAC);
            res.put("GATS" + lag, GearyAC);
        }

        return res;
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

    public boolean isLtype(int atomIdx, IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        IAtom CurAt =  m.getAtom(atomIdx);
        boolean tL=false;

        // Hydrogens
        int H = 0;
        try {
            H = CurAt.getImplicitHydrogenCount();
        } catch (Exception e) { }

        // counters
        int VD = 0;
        int nC = 0, nNonC = 0;
        for (int j=0; j<nSK; j++) {
            if (j==atomIdx) continue;
            if (ConnAugMatrix[atomIdx][j]>0) {

                VD++;

                if (ConnAugMatrix[j][j] == 6)
                    nC++;
                else if (ConnAugMatrix[j][j] != 1)
                    nNonC++;


            }

        }

        if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
            tL = true;

        if (CurAt.getSymbol().equalsIgnoreCase("Br"))
            tL = true;

        if (CurAt.getSymbol().equalsIgnoreCase("I"))
            tL = true;

        if (CurAt.getSymbol().equalsIgnoreCase("C")) {
            if ( (nC > 0) && (nNonC == 0) )
                tL = true;
        }

        if (CurAt.getSymbol().equalsIgnoreCase("S")) {
            if ( (nC == 2) && ( (VD+H) == 2 ))
                tL = true;
        }

        return tL;
    }

}
