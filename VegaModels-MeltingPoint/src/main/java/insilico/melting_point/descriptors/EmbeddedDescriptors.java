package insilico.melting_point.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.Rings;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.other.WeightsHydrophobicityGC;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.melting_point.descriptors.weights.MoleculePaths;
//import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.*;

import java.util.*;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;

    public double piPC01 = MISSING_VALUE;
    public double CATS2D_00_DD = MISSING_VALUE;
    public double ATS2i = MISSING_VALUE;
    public double P_VSA_LogP_3 = MISSING_VALUE;
    public double nCIC = MISSING_VALUE;
    public double CIC3 = MISSING_VALUE;
    public double GATS1e = MISSING_VALUE;
    public double BIC1 = MISSING_VALUE;
    public double ATSC1i = MISSING_VALUE;
    public double nArCOOH = MISSING_VALUE;
    public double NPerc = MISSING_VALUE;
    public double B02_C_O= MISSING_VALUE;
    public double ATS3m= MISSING_VALUE;
    public double nCconj= MISSING_VALUE;
    public double F04_C_O= MISSING_VALUE;
    public double nROH= MISSING_VALUE;
    public double CATS2D_04_PL = MISSING_VALUE;
    public double S_106= MISSING_VALUE;
    public double T_N_N= MISSING_VALUE;
    public double N_072= MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws GenericFailureException {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws GenericFailureException {

        CalculatePiPC(mol);
        CalculateCATS2DDD(mol);
        CalculateCATS2DPL(mol);
        CalculateAutoCorrelation(mol);
        CalculatePVSA(mol);
        CalculateNCic(mol);
        CalculateIC(mol);
        CalculateFunctionalGroups(mol);
        CalculateConstitutional(mol);
        CalculateACF(mol);
        CalculateAtomPairs2D(mol);
    }

    private void CalculateAtomPairs2D(InsilicoMolecule mol) {

        String[][] ATOM_COUPLES = {
                {"N", "N"},
                {"C", "O"}
        };

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
            log.warn("Invalid structure, unable to calculate: " + mol.getInputSMILES());
            return;
        }

        for (int d = 0; d< ATOM_COUPLES.length; d++) {

            int descT = 0;
            int[] descB = new int[10];
            int[] descF = new int[10];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 1) // Consistent with D7 - ignore atoms directly attached
                                descT += TopoMat[i][j];

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= 10) {
                                descB[TopoMat[i][j]-1] = 1;
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (ATOM_COUPLES[d][0].compareTo(ATOM_COUPLES[d][1]) == 0) {
                descT /= 2;
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            if(Arrays.equals(ATOM_COUPLES[d], new String[]{"N", "N"}))
                T_N_N = descT;

            for (int i=0; i<descB.length; i++) {
                int lag = i+1;
                if(lag == 2 && Arrays.equals(ATOM_COUPLES[d], new String[]{"C", "O"}))
                    B02_C_O = descB[i];
                if(lag == 4 && Arrays.equals(ATOM_COUPLES[d], new String[]{"C", "O"}))
                    F04_C_O = descF[i];

            }
        }
    }

    private void CalculateACF(InsilicoMolecule mol) {
        try {
            AtomCenteredFragments block = new AtomCenteredFragments();
            block.Calculate(mol);
            N_072 = block.GetByName("N-072").getValue();
            S_106 = block.GetByName("S-106").getValue();
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }

    private void CalculateConstitutional(InsilicoMolecule mol) {
        try {
            Constitutional block = new Constitutional();
            block.Calculate(mol);
            NPerc = block.GetByName("NPerc").getValue();
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }

    }

    private void CalculateFunctionalGroups(InsilicoMolecule mol) {

        try {
            FunctionalGroups block = new FunctionalGroups();
            block.Calculate(mol);
            nArCOOH = block.GetByName("nArCOOH").getValue();
            nCconj = block.GetByName("nCconj").getValue();
            nROH = block.GetByName("nROH").getValue();
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }

    private void CalculateIC(InsilicoMolecule mol) {


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
            ConnMatH = ConnectionAugMatrix.getMatrix(curMolH);
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSKH = curMolH.getAtomCount();
        int[] VertexDegH = WeightsVertexDegree.getWeights(curMolH, false);

        double bic_denom = 0;
        for (IBond bnd : curMolH.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag = 1; CurLag <= 8; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSKH);
            for (int i = 0; i < nSKH; i++) {
                IAtom atStart = curMolH.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j = 0; j < nSKH; j++) {
                    if (i == j) continue;
                    if (TopoDistMatH[i][j] == CurLag) {
                        IAtom atEnd = curMolH.getAtom(j);
                        ShortestPaths shortestPaths = new ShortestPaths(curMolH, atStart);
                        List<IAtom> sp = Arrays.asList(shortestPaths.atomsTo(atEnd));
                        StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                        for (int k = 0; k < (sp.size() - 1); k++) {
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
            for (int i = 0; i < nSKH; i++) {
                ArrayList<String> CurNeig = NeigList.get(i);
                boolean foundMatch = false;
                for (int k = 0; k < G.size(); k++) {
                    if (CompareNeigVector(CurNeig, G.get(k))) {
                        foundMatch = true;
                        int buf = Gn.get(k);
                        Gn.set(k, (buf + 1));
                        break;
                    }
                }
                if (!foundMatch) {
                    G.add(CurNeig);
                    Gn.add(1);
                }
            }

            // Calculate IC and CIC indices

            double ic = 0;
            for (int i = 0; i < Gn.size(); i++)
                ic += ((double) Gn.get(i) / nSKH) * (Math.log((double) Gn.get(i) / nSKH));
            ic = (-1.00 / Math.log(2)) * ic;

            double diff = Math.log(nSKH) / Math.log(2);

            if(CurLag == 1)
                BIC1 = bic_denom == 0 ? 0 : ic / bic_denom;

            if(CurLag == 3)
                CIC3 =  diff - ic;
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


    private void CalculateNCic(InsilicoMolecule mol) {

        try {
            Rings block = new Rings();
            block.Calculate(mol);
            nCIC = block.GetByName("nCIC").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

        try {

            // get ring sets directly from molecule (cache)
            IRingSet sssr = mol.GetSSSR();
            nCIC = sssr.getAtomContainerCount();


        } catch (Throwable e) {
            log.warn(e.getMessage());
        }

    }

    private void CalculatePVSA(InsilicoMolecule mol) {


        IAtomContainer m = null;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.getInputSMILES());
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
            double refR = GetRefBondLength(bo.getAtom(0), bo.getAtom(1));
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

        double[] w = (new WeightsHydrophobicityGC()).getWeightsForFragmentId(ACF);

        int b = 3;
        double PVSA = 0;
        for (int i = 0; i < nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE) continue;
            if (b == CalculateBin(w[i]))
                PVSA += VSA[i];
        }

        P_VSA_LogP_3 = PVSA;

    }

    private void CalculateAutoCorrelation(InsilicoMolecule mol) {
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.getInputSMILES());
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

        // pesi i
        double[] w = (new WeightsIonizationPotential()).getScaledWeights(m);
        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        for (int lag=1; lag<=10; lag++) {

            if(lag > 2)
                break;

            double AC=0, ACS=0;

            for (int i=0; i<nSK; i++) {

//                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        AC += w[i] * w[j];
                        ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                    }
            }

            // AC transformed in log form
            AC /= 2.0;
            AC = Math.log(1 + AC);

            ACS /= 2.0;

            if(lag == 1)
                ATSC1i = ACS;
            if(lag == 2)
                ATS2i = AC;
        }


        w = (new WeightsElectronegativity()).getScaledWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        for (int lag=1; lag<=10; lag++) {

            if(lag > 1)
                break;

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


            GATS1e = GearyAC;
        }


        w = (new WeightsMass()).getScaledWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        for (int lag=1; lag<=10; lag++) {

            if(lag < 3)
                continue;
            if(lag > 3)
                break;

            double AC=0;

            for (int i=0; i<nSK; i++) {

//                denom += Math.pow((w[i] - wA), 2);

                for (int j=0; j<nSK; j++)
                    if (TopoMatrix[i][j] == lag) {
                        AC += w[i] * w[j];
                    }
            }



            // AC transformed in log form
            AC /= 2.0;
            AC = Math.log(1 + AC);

            ATS3m = AC;
        }



    }

    private void CalculateCATS2DPL(InsilicoMolecule mol){

        String[][][] AtomCouples = {{{"P"}, {"L"}}};


        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.getInputSMILES());
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


        for (String[][] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
//                if (isIn(CatsTypes[i], atomCouple[0][0])) {
                for (int j = i; j < nSK; j++) {
//                        if (isIn(CatsTypes[j], atomCouple[1][0])) {

                    if (isInCatsType(CatsTypes[i], CatsTypes[j], atomCouple[0][0], atomCouple[1][0]))
                        if (TopoMat[i][j] < 10)
                            desc[TopoMat[i][j]]++;

//                        }
                }
//                }
            }

            for (int i = 0; i < desc.length; i++) {
                if(i == 4)
                    CATS2D_04_PL = desc[i];
                else continue;
            }
        }

    }

    private void CalculateCATS2DDD(InsilicoMolecule mol) {


        String[][][] AtomCouples = {{{"D"}, {"D"}}};

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure: " + mol.getInputSMILES());
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

        for (String[][] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
//                if (isIn(CatsTypes[i], atomCouple[0][0])) {
                for (int j = i; j < nSK; j++) {
//                        if (isIn(CatsTypes[j], atomCouple[1][0])) {

                    if (isInCatsType(CatsTypes[i], CatsTypes[j], atomCouple[0][0], atomCouple[1][0]))
                        if (TopoMat[i][j] < 10)
                            desc[TopoMat[i][j]]++;

//                        }
                }
//                }
            }

            for (int i = 0; i < desc.length; i++) {
                if(i == 0)
                    CATS2D_00_DD = desc[i];
            }
        }
    }

    private void CalculatePiPC(InsilicoMolecule Mol) throws GenericFailureException {
        MoleculePaths paths = new MoleculePaths(Mol);
        piPC01 = Math.log(1 + paths.Multiple_Path_Counts[0]);
    }

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

    private int CalculateBin(double value) {

        // LOGP
        if (value < -1.5) return 1;
        if (value < -0.5) return 2;
        if (value < -0.25) return 3;
        if (value < 0) return 4;
        if (value < 0.25) return 5;
        if (value < 0.52) return 6;
        if (value < 0.75) return 7;
        return 8;
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
    private

    boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
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

        return -999;
    }


}
