package insilico.tpo_oberon.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.weights.basic.*;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsQuantumNumber;
import insilico.core.descriptor.blocks.weights.other.WeightsValenceVertexDegree;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.tpo_oberon.descriptors.weights.DescriptorMLogP;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;

import java.util.ArrayList;
import java.util.Arrays;


public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private final double MISSING_VALUE = -999;

    public double GATS1e = MISSING_VALUE;
    public double nArOH = MISSING_VALUE;
    public double CATS2D_02_DL = MISSING_VALUE;
    public double MATS1e = MISSING_VALUE;
    public double MATS1s = MISSING_VALUE;
    public double C_026 = MISSING_VALUE;
    public double CATS2D_03_DL = MISSING_VALUE;
    public double B10_C_C = MISSING_VALUE;
    public double MATS1p = MISSING_VALUE;
    public double nCb_ = MISSING_VALUE;
    public double nX = MISSING_VALUE;
    public double Uc = MISSING_VALUE;
    public double P_VSA_i_1 = MISSING_VALUE;
    public double SpMAD_B_v_ = MISSING_VALUE;
    public double nCbH = MISSING_VALUE;
    public double GATS1s = MISSING_VALUE;
    public double MATS1m = MISSING_VALUE;
    public double MLOGP = MISSING_VALUE;
    public double SpMax2_Bh_s = MISSING_VALUE;
    public double Eta_C_A = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws DescriptorNotFoundException, GenericFailureException {
        CalculateDescriptors(mol);

    }

    private void CalculateDescriptors(InsilicoMolecule mol) throws DescriptorNotFoundException, GenericFailureException {

        CalculateAutocorrelation(mol);
        CalculateFunctionalGroups(mol);
        CalculateCATS2D(mol);
        CalculateACF(mol);
        CalculateAtomPairs2D(mol);
        CalculateConstitutional(mol);
        CalculateMoleculearProperties(mol);
        CalculatePVSA(mol);
        CalculateMatrices2D(mol);
        CalculateBurdenEigenValues(mol);
        CalculateMLOGP(mol);
        CalculateEtaIndeces(mol);
    }

    private void CalculateEtaIndeces(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            Eta_C_A = MISSING_VALUE;
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
            Eta_C_A = MISSING_VALUE;
            return;
        }

        RingSet MolRings;
        try {
            MolRings = mol.GetSSSR();
        } catch (Exception e) {
            log.warn(e.getMessage());
            Eta_C_A = MISSING_VALUE;
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
                Eta_C_A = MISSING_VALUE;
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
            gamma[i] = alpha[i] / (beta_s[i] + beta_ns[i]);
            beta_s[i] /= 2.0;
            beta_ns[i] = (beta_ns[i] - correctionFactor) / 2.0 + correctionFactor;

        }

        //// ETA indices
        // composite index
        double eta_c = 0;
        for (int i=0; i<(nSK-1); i++)
            for (int j=(i+1); j<nSK; j++)
                eta_c += Math.pow( (gamma[i] * gamma[j]) / Math.pow(TopoMatrix[i][j],2) , 0.5);

        Eta_C_A = eta_c / (double)nSK;
    }

    private void CalculateMLOGP(InsilicoMolecule mol) {

        DescriptorMLogP descriptorMLogP = new DescriptorMLogP();
        MLOGP = descriptorMLogP.Calculate(mol);
    }

    private void CalculateBurdenEigenValues(InsilicoMolecule mol) {


        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure " + mol.GetSMILES());
            return;
        }

        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();


        WeightsIState weightsIState = new WeightsIState();
        double[] w = weightsIState.getWeights(m, true);

        // correction for compatibility with D7
        // H I-state is always 1
        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }

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

        for (int i=0; i<8; i++) {
            if(i == 1){
                double valH;
                if (i>(eigenvalues.length-1)) {
                    valH = 0;
                } else {
                    if (eigenvalues[eigenvalues.length-1-i] > 0)
                        valH = eigenvalues[eigenvalues.length-1-i];
                    else
                        valH = 0;
                }
                SpMax2_Bh_s = valH;
            } else { }
        }

    }

    private void CalculateMatrices2D(InsilicoMolecule mol) throws GenericFailureException {


        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid Structure: " + mol.GetSMILES());
            return;
        }

        // Adj matrix available for calculations

        int nSK = curMol.getAtomCount();

        double[][] Mat = mol.GetMatrixBurden();
        double[] w = ( new WeightsVanDerWaals()).getScaledWeights(curMol);

        for (int i=0; i<nSK; i++)
            Mat[i][i] = w[i];

        Matrix DataMatrix = new Matrix(Mat);
        double[] eigenvalues;

        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            log.warn("Unable to calculate eigenvalues: " + e.getMessage());
            return;
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
        EigAve = EigAve / (double) nSK;
        double EigDev = 0;
        for (double val : eigenvalues)
            EigDev += Math.abs(val - EigAve);
        double NormEigDev = EigDev / (double) nSK;



        // Eigenvalue spectral moments
        double[] SpecMoments = new double[7];
        for (double d : SpecMoments) d = 0;
        for (double val : eigenvalues) {
            for (int expIdx=1; expIdx<7; expIdx++)
                SpecMoments[expIdx] += Math.pow(val, expIdx);
        }
        SpecMoments[1] = Math.signum(SpecMoments[1]) * Math.log(1 + Math.abs(SpecMoments[1]));
        SpecMoments[2] = Math.log(1 + SpecMoments[2]);
        SpecMoments[3] = Math.signum(SpecMoments[3]) * Math.log(1 + Math.abs(SpecMoments[3]));
        SpecMoments[4] = Math.log(1 + SpecMoments[4]);
        SpecMoments[5] = Math.signum(SpecMoments[5]) * Math.log(1 + Math.abs(SpecMoments[5]));
        SpecMoments[6] = Math.log(1 + SpecMoments[6]);

        SpMAD_B_v_ = NormEigDev;


    }

    private void CalculatePVSA(InsilicoMolecule mol) {


        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid strucutre " + mol.GetSMILES());
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

        double[] w = (new WeightsIonizationPotential()).getScaledWeights(m);
        for (int b=1; b<=4; b++) {

            if(b == 1) {
                double PVSA = 0;


                for (int i = 0; i < nSK; i++) {
                    if (w[i] == Descriptor.MISSING_VALUE) continue;
                    if (b == CalculateBin(w[i]))
                        PVSA += VSA[i];
                }


                P_VSA_i_1 = PVSA;
            }
        }
    }

    private void CalculateMoleculearProperties(InsilicoMolecule mol) {

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure " + mol.GetSMILES());
            return;
        }

        int nSK = curMol.getAtomCount();
        int nBO = curMol.getBondCount();

        // Gets matrices
        double[][] AdjConnMatrix;
        try {
            AdjConnMatrix = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

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
        // Correction for NO2 - Dragon 7 compliance
        for (int i=0; i<nSK; i++) {

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
                        nDblBonds++;
                        scbo += (+2 - 1);

                    }
                }
            }
        }

        Uc = Math.log(1 + nDblBonds + nTrpBonds + nArBonds) / Math.log(2);
    }

    private void CalculateConstitutional(InsilicoMolecule mol) {

        try {
            DescriptorBlock block = new Constitutional();
            block.Calculate(mol);
            nX = block.GetByName("nX").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateAtomPairs2D(InsilicoMolecule mol) {

        String[][] ATOM_COUPLES = {{"C", "C"}};

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
            log.warn("Invalid structure, unable to calculate: " + mol.GetSMILES());
            return;
        }

        for (int d = 0; d< ATOM_COUPLES.length; d++) {

            int descT = 0;
            int[] descB = new int[10];
            int[] descF = new int[10];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);

            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2) // DA VEDERE PERCHE MAGGIORE DI 2
                                descT += TopoMat[i][j];

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= 10) {
                                descB[TopoMat[i][j] - 1] = 1;
                                descF[TopoMat[i][j] - 1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (ATOM_COUPLES[d][0].compareTo(ATOM_COUPLES[d][1]) == 0) {
                descT /= 2;
                for (int i = 0; i < descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i=0; i<descB.length; i++) {
                int lag = i + 1;
                if(lag == 10)
                    B10_C_C = descB[i];
            }


        }
    }

    private void CalculateACF(InsilicoMolecule mol) {
        try {
            DescriptorBlock block = new AtomCenteredFragments();
            block.Calculate(mol);
            C_026 = block.GetByName("C-026").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateCATS2D(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure " + mol.GetSMILES());
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

                }
            }

            for (int i = 0; i < desc.length; i++) {
                if(i==2)
                    CATS2D_02_DL = desc[i];
                if(i==3)
                    CATS2D_03_DL = desc[i];
            }
        }


    }

    private void CalculateFunctionalGroups(InsilicoMolecule mol) {

        try {
            DescriptorBlock block = new FunctionalGroups();
            block.Calculate(mol);
            nArOH = block.GetByName("nArOH").getValue();
            nCb_ = block.GetByName("nCb–").getValue();
            nCbH = block.GetByName("nCbH").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateAutocorrelation(InsilicoMolecule mol) {
        CalculateGATS(mol);
        CalculateMATS(mol);
    }

    private void CalculateMATS(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
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

        // s
        double[] w = (new WeightsIState()).getWeights(m, true);
        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        for (int lag=1; lag<=8; lag++) {

            if(lag == 1){
                double MoranAC=0;
                double denom = 0, delta = 0;

                for (int i=0; i<nSK; i++) {

                    denom += Math.pow((w[i] - wA), 2);

                    for (int j=0; j<nSK; j++)
                        if (TopoMatrix[i][j] == lag) {
                            MoranAC += (w[i] - wA) * (w[j] - wA);
                            delta++;
                        }
                }

                if (delta > 0) {
                    if (denom == 0) {
                        MATS1s = 1;
                    } else {
                        MATS1s = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                    }
                }
            }
        }

        // s
        w = (new WeightsElectronegativity()).getWeights(m);

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        for (int lag=1; lag<=8; lag++) {

            if(lag == 1){
                double MoranAC=0;
                double denom = 0, delta = 0;

                for (int i=0; i<nSK; i++) {

                    denom += Math.pow((w[i] - wA), 2);

                    for (int j=0; j<nSK; j++)
                        if (TopoMatrix[i][j] == lag) {
                            MoranAC += (w[i] - wA) * (w[j] - wA);
                            delta++;
                        }
                }

                if (delta > 0) {
                    if (denom == 0) {
                        MATS1e = 1;
                    } else {
                        MATS1e = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                    }
                }
            }
        }

        // m
        w = (new WeightsMass()).getWeights(m);

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        for (int lag=1; lag<=8; lag++) {

            if(lag == 1){
                double MoranAC=0;
                double denom = 0, delta = 0;

                for (int i=0; i<nSK; i++) {

                    denom += Math.pow((w[i] - wA), 2);

                    for (int j=0; j<nSK; j++)
                        if (TopoMatrix[i][j] == lag) {
                            MoranAC += (w[i] - wA) * (w[j] - wA);
                            delta++;
                        }
                }

                if (delta > 0) {
                    if (denom == 0) {
                        MATS1m = 1;
                    } else {
                        MATS1m = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                    }
                }
            }
        }

        // p
        w = (new WeightsPolarizability()).getWeights(m);

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        for (int lag=1; lag<=8; lag++) {

            if(lag == 1){
                double MoranAC=0;
                double denom = 0, delta = 0;

                for (int i=0; i<nSK; i++) {

                    denom += Math.pow((w[i] - wA), 2);

                    for (int j=0; j<nSK; j++)
                        if (TopoMatrix[i][j] == lag) {
                            MoranAC += (w[i] - wA) * (w[j] - wA);
                            delta++;
                        }
                }

                if (delta > 0) {
                    if (denom == 0) {
                        MATS1p = 1;
                    } else {
                        MATS1p = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                    }
                }
            }
        }

    }

    private void CalculateGATS(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn("Invalid structure: " + mol.GetSMILES());
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

        // s
        double[] w = (new WeightsIState()).getWeights(m, true);
        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }

        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        // Calculates autocorrelations
        for (int lag=1; lag<=8; lag++) {

            if(lag==1){
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
                    if (denom !=0){
                        GATS1s = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double) (nSK - 1))) * denom);
                    }
                }
            }
        }


        w = (new WeightsElectronegativity()).getWeights(m);

        // Calculates weights averages
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        // Calculates autocorrelations
        for (int lag=1; lag<=8; lag++) {

            if(lag==1){
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
                    if (denom !=0){
                        GATS1e = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double) (nSK - 1))) * denom);
                    }
                }
            }
        }
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

    private int CalculateBin(double value) {

        if (value < 1.0) return 1;
        if (value < 1.5) return 2;
        if (value < 2.0) return 3;
        if (value < 2.5) return 4;
        if (value < 3.0) return 5;
        return 6;
    }

    private final String[] TYPE_D = { "D", ""} ;
    private final String[] TYPE_A = { "A", ""};
    private final String[] TYPE_P = { "P", ""};
    private final String[] TYPE_N = { "N", ""};
    private final String[] TYPE_L = { "L", ""};
    private final String[] TYPE_CYC = { "Cyc", ""};

    private final String[][][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };

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



            if (tD) AtomTypes[i].add(TYPE_D[0]);
            if (tL) AtomTypes[i].add(TYPE_L[0]);

        }

        return AtomTypes;
    }


}
