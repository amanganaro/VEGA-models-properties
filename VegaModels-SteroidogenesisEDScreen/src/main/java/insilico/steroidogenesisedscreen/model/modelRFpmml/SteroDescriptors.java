package insilico.steroidogenesisedscreen.model.modelRFpmml;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.descriptor.blocks.weights.other.WeightsQuantumNumber;
import insilico.core.descriptor.blocks.weights.other.WeightsValenceVertexDegree;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;

import java.util.*;

public class SteroDescriptors {

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
        double[][] DistDetMatrix;
        double[][] BurdenMat_s;
        double[][] AdjConnMatrix;
        double[][][] EdgeAdjMat;
        try {
            TopoMatrix =  mol.GetMatrixTopologicalDistance();
            DistDetMatrix = mol.GetMatrixDistanceDetour();
            BurdenMat_s = BurdenMatrix.getMatrix(m);
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
            AdjConnMatrix = mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate topological distance matrix");
        }

        double[] w = (new WeightsIState()).getWeights(m, false);
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing a E-state weight value");
        // Builds the weighted matrix
        for (int i=0; i<nSK; i++)
            BurdenMat_s[i][i] = w[i];

        RingSet MolRings;
        try {
            MolRings = mol.GetSSSR();
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate ring structures");
        }


        //// MAXDP & SaaCH

        double D01_MaxDP = Descriptor.MISSING_VALUE;


        // Get I-States weights
        WeightsIState w_istate = new WeightsIState();
        double[] w_is = w_istate.getWeights(m, false);
        for (double val : w_is)
            if (val == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("unable to calculate all E-states");

        // Calculate E-States
        double[] w_es = new double[nSK];
        for (int at = 0; at<nSK; at++) {
            double sumDeltaI = 0;

            for (int j = 0; j < nSK; j++)
                if (at != j)
                    sumDeltaI += (w_is[at] - w_is[j]) / Math.pow((double) TopoMatrix[at][j] + 1.0, 2.0);

            w_es[at] = w_is[at] + sumDeltaI;

            // max variation
            if (D01_MaxDP == Descriptor.MISSING_VALUE)
                D01_MaxDP = sumDeltaI;
            if (sumDeltaI > D01_MaxDP)
                D01_MaxDP = sumDeltaI;
        }

        double D17_SaaCH = 0;
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

            // aaCH group
            if ( (nBnd==2) && (nAr == 2) && (nH == 1) && (Charge == 0) )
                D17_SaaCH += w_es[at];
        }


        //// IDE

        // Information content
        int[] TopoDistFreq = new int[nSK];  // frequencies of topological distances
        double TopoDistFreqSum = 0;
        for (int i=0; i<nSK; i++)
            TopoDistFreq[i] = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i+1; j<nSK; j++) {
                TopoDistFreq[TopoMatrix[i][j]]++;
                TopoDistFreqSum += TopoMatrix[i][j];
            }

        // IDE
        double D02_IDE = 0;
        double denom = (double)nSK*(nSK-1)/2.00;
        for (int i=0; i<nSK; i++)
            if (TopoDistFreq[i]>0)
                D02_IDE += ((double)TopoDistFreq[i]/denom) * (Math.log((double)TopoDistFreq[i]/denom));
        D02_IDE = (-1.00 / Math.log(2)) * D02_IDE;


        //// SpMaxA_D/Dt

        double D04_SpMaxA_D_Dt = CalculateEigDescriptors(DistDetMatrix, nSK).get("SpMaxA");


        //// SpMax_D(s)

        double D06_SpMax_B_s = CalculateEigDescriptors(BurdenMat_s, nSK).get("SpMax");





        //// Edge adjacency

        double[][] ea_mat = new double[nBO][nBO];
        w = new double[nBO];

        // EA(bo)
        for (int i=0; i<nBO; i++) {
            w[i] = MoleculeUtilities.Bond2Double(m.getBond(i));
            for (int j = 0; j < nBO; j++)
                ea_mat[i][j] = EdgeAdjMat[i][j][0];
        }
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                if (ea_mat[i][j] != 0)
                    ea_mat[i][j] = w[j];
        double D13_SM05_EA_bo = CalculateEigDescriptors(ea_mat, nBO).get("SM05");


        // EA(dm)
        for (int i=0; i<nBO; i++)
            for (int j = 0; j < nBO; j++)
                ea_mat[i][j] = EdgeAdjMat[i][j][0];
        for (int i=0; i<nBO; i++) {
            IAtom a =  m.getBond(i).getAtom(0);
            IAtom b =  m.getBond(i).getAtom(1);
            double CurVal = GetDipoleMoment(m, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(m, b, a);
            w[i] = CurVal;
        }
        for (int i=0; i<nBO; i++)
            for (int j=0; j<nBO; j++)
                if (ea_mat[i][j] != 0)
                    ea_mat[i][j] = w[j];
        double D14_SM02_EA_dm = CalculateEigDescriptors(ea_mat, nBO).get("SM02");


        // AEA(bo)
        for (int i=0; i<nBO; i++) {
            w[i] = MoleculeUtilities.Bond2Double(m.getBond(i));
            for (int j = 0; j < nBO; j++)
                ea_mat[i][j] = EdgeAdjMat[i][j][0];
        }
        for (int i=0; i<nBO; i++)
            ea_mat[i][i] = w[i];
        double D16_EIG12_AEA_bo = CalculateEigDescriptors(ea_mat, nBO).get("EIG12");


        // ETA

        int[] AtomicNumber = new int[nSK];
        int[] ElectronNumber = new int[nSK]; // valence electron number
        int[] QuantumNumber = new int[nSK];
        int[] VertexDegree = WeightsVertexDegree.getWeights(m, true); // include H in VD
        boolean[] isAromatic = new boolean[nSK];

        WeightsValenceVertexDegree wVVD = new WeightsValenceVertexDegree();
        WeightsQuantumNumber wQN = new WeightsQuantumNumber();
        for (int i=0; i<nSK; i++) {
            AtomicNumber[i] = m.getAtom(i).getAtomicNumber();
            ElectronNumber[i] = wVVD.GetValenceElectronsNumber(m.getAtom(i).getSymbol());
            QuantumNumber[i] = wQN.getWeight(m.getAtom(i).getSymbol());
            if ( (ElectronNumber[i] == Descriptor.MISSING_VALUE) || (QuantumNumber[i] == Descriptor.MISSING_VALUE) )
                throw new GenericFailureException("unable to calculate electron or quantum number");
            isAromatic[i] = m.getAtom(i).getFlag(CDKConstants.ISAROMATIC);
        }

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
                boolean exception = m.getAtom(i).getSymbol().equalsIgnoreCase("N") && (VertexDegree[i] == 3);
                // P in pyrrole
                exception = exception || (m.getAtom(i).getSymbol().equalsIgnoreCase("P") && (VertexDegree[i] == 3));
                // S, O, Se
                exception = exception || (m.getAtom(i).getSymbol().equalsIgnoreCase("O") ||
                        m.getAtom(i).getSymbol().equalsIgnoreCase("S") || m.getAtom(i).getSymbol().equalsIgnoreCase("Se") );

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
                    TotBondOrder += Manipulator.CountImplicitHydrogens(m.getAtom(i));

                    if ( (m.getAtom(i).getSymbol().equalsIgnoreCase("O") ) || (m.getAtom(i).getSymbol().equalsIgnoreCase("S") ) )
                        if (TotBondOrder == 2)
                            correctionFactor = 0.5;
                    if ( (m.getAtom(i).getSymbol().equalsIgnoreCase("N") ) || (m.getAtom(i).getSymbol().equalsIgnoreCase("P") ) )
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
                        if ( (curRing.contains(m.getAtom(i))) && (curRing.contains(m.getAtom(at_idx))) ) {
                            if ( (m.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) && (m.getAtom(at_idx).getFlag(CDKConstants.ISAROMATIC)) ) {
                                SameAromRing = true;
                                break;
                            }
                        }
                    }

                    if ( (!SameAromRing) && (isConjugated(m, i, AdjConnMatrix)) && (isConjugated(m, at_idx, AdjConnMatrix)) ) {

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
                if (m.getAtom(i).getFormalCharge() == +1) {
                    int Odbl = 0, Osngminus = 0, vd = 0;
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (AdjConnMatrix[i][j] != 0) {
                            vd++;
                            if (AdjConnMatrix[j][j] == 8) {
                                if (AdjConnMatrix[i][j] == 2)
                                    Odbl++;
                                if (AdjConnMatrix[i][j] == 1)
                                    if (m.getAtom(j).getFormalCharge() == -1)
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

        double eta_betas = 0;
        double eta_betans = 0;
        for (int i=0; i<nSK; i++) {
            eta_betas += beta_s[i];
            eta_betans += beta_ns[i];
        }

        double eta_beta = eta_betas + eta_betans;
        double D12_ETA_BETA_A = eta_beta / (double) nSK;




        //////// descriptors with H-filled structure

        IAtomContainer curMolH = Manipulator.AddHydrogens(m);

        double[][] ConnMatH;
        int[][] TopoDistMatH;
        double[][] BurdenMat_h;
        double[][] AdjConnMatrix_h;
        try {
            ConnMatH= ConnectionAugMatrix.getMatrix(curMolH);
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
            BurdenMat_h = BurdenMatrix.getMatrix(curMolH);
            AdjConnMatrix_h = ConnectionAugMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate matrix - " + e);
        }
        int nSKH = curMolH.getAtomCount();
        int nBOH = curMolH.getBondCount();


        //// SIC1

        int[] VertexDegH = WeightsVertexDegree.getWeights(curMolH, false);
        int CurLag=1;

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

        double ic=0;
        for (int i=0; i<Gn.size(); i++)
            ic += ((double)Gn.get(i)/nSKH) * (Math.log((double)Gn.get(i)/nSKH));
        ic = (-1.00 / Math.log(2)) * ic;

        double diff = Math.log(nSKH) / Math.log(2);
        double D03_SIC1 = ic / diff;


        //// SpMax2_Bh(s)

        w = (new WeightsIState()).getWeights(curMolH, true);
        // correction for compatibility with D7
        // H I-state is always 1
        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing a E-state weight value");
        // Builds the weighted matrix
        for (int i=0; i<nSK; i++)
            BurdenMat_h[i][i] = w[i];

        double D08_SpMax2_Bh_s = CalculateEigDescriptors(BurdenMat_h, nSKH).get("SpMax2");


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

        double D10_P_VSA_ppp_ar = 0;
        double D11_P_VSA_ppp_con = 0;
        for (int i = 0; i < nSK; i++) {
            if (curMolH.getAtom(i).isAromatic())
                D10_P_VSA_ppp_ar += VSA[i];
            if (isConjugated(curMolH, i, AdjConnMatrix_h))
                D11_P_VSA_ppp_con += VSA[i];
        }


        // GATS

        w = (new WeightsIonizationPotential()).getScaledWeights(curMolH);

        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing ionization potential weight");

        double wA = 0;
        for (int i=0; i<nSKH; i++)
            wA += w[i];
        wA = wA / ((double) nSKH);

        int lag = 7;
        double D07_GATS7i = 0;
        double ac_denom = 0, delta = 0;

        for (int i=0; i<nSKH; i++) {
            ac_denom += Math.pow((w[i] - wA), 2);
            for (int j=0; j<nSKH; j++)
                if (TopoDistMatH[i][j] == lag) {
                    D07_GATS7i += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (ac_denom == 0)
                D07_GATS7i = 0;
            else
                D07_GATS7i = ((1 / (2 * delta)) * D07_GATS7i) / ((1 / ((double)(nSKH - 1))) * ac_denom);
        }



        double[] res = new double[14];
        res[0] = D01_MaxDP;
        res[1] = D02_IDE;
        res[2] = D03_SIC1;
        res[3] = D04_SpMaxA_D_Dt;
        res[4] = D06_SpMax_B_s;
        res[5] = D07_GATS7i;
        res[6] = D08_SpMax2_Bh_s;
        res[7] = D10_P_VSA_ppp_ar;
        res[8] = D11_P_VSA_ppp_con;
        res[9] = D12_ETA_BETA_A;
        res[10] = D13_SM05_EA_bo;
        res[11] = D14_SM02_EA_dm;
        res[12] = D16_EIG12_AEA_bo;
        res[13] = D17_SaaCH;

        return res;

    }

    public final static String[] DESCRIPTOR_NAMES = {
        "MAXDP",
        "IDE",
        "SIC1",
        "SpMaxA_D/Dt",
        "SpMax_B(s)",
        "GATS7i",
        "SpMax2_Bh(s)",
        "P_VSA_ppp_ar",
        "P_VSA_ppp_con",
        "Eta_beta_A",
        "SM05_EA(bo)",
        "SM02_EA(dm)",
        "Eig12_AEA(bo)",
        "SaaCH",
    };

    // values calculated from the saved TRAINING SET only
    private final double[] DESC_AVE =
            {3.562355, 2.797488, 0.573561, 0.674663, 6.783517, 0.690909, 5.809250, 39.761580, 51.654524, 1.046637, 6.347451, 1.761847, 0.303109, 5.257276};
    private final double[] DESC_STD =
            {1.556581, 0.662891, 0.140600, 0.177682, 1.190133, 0.485572, 1.344640, 42.510800, 50.327861, 0.304315, 1.781037, 1.222392, 0.785889, 6.268043};

    public double[] Scale(double[] descriptors) {
        double[] scaled = new double[descriptors.length];
        for (int i=0; i<descriptors.length; i++)
            scaled[i] = ( descriptors[i] - DESC_AVE[i] ) / DESC_STD[i];
        return scaled;
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

    private HashMap<String, Double> CalculateEigDescriptors(double[][] mat, int nSK) throws GenericFailureException {
        HashMap<String, Double> res = new HashMap<>();

        Matrix DataMatrix = new Matrix(mat);
        double[] eigenvalues;
        Matrix eigenvectors;
        try {
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            eigenvectors = ed.getV();
            Arrays.sort(eigenvalues);
        } catch (Throwable e) {
            throw new GenericFailureException("unable to calculate eigenvalues - " + e);
        }

        // check on precision
        for (int i=0; i<eigenvalues.length; i++) {
            if (Math.abs(eigenvalues[i]) < 0.000001)
                eigenvalues[i] = 0;
        }

        double EigMax = eigenvalues[0];
        for (double val : eigenvalues)
            if (val > EigMax)
                EigMax = val;

        res.put("SpMax", EigMax);
        res.put("SpMaxA", EigMax / (double) nSK);
        double SpMax2 = nSK>=2 ? eigenvalues[eigenvalues.length - 2] : 0;
        SpMax2 = SpMax2 > 0 ? SpMax2 : 0;
        res.put("SpMax2", SpMax2);


        // Eig
        int idx_12 = eigenvalues.length - 12;
        double eig12 = idx_12 >=0 ? eigenvalues[idx_12] : 0;
        res.put("EIG12", eig12);

        int idx_14 = eigenvalues.length - 14;
        double eig14 = idx_14 >=0 ? eigenvalues[idx_14] : 0;
        res.put("EIG14", eig14);


        // SM
        double[] SpecMoments = new double[8];
        for (int i=0; i<8; i++) {
            SpecMoments[i] = 0;
            for (double val : eigenvalues) {
                if (Math.abs(val) > 0)
                    SpecMoments[i] += Math.pow(val, (i + 1));
            }
            SpecMoments[i] = Math.log(1 + SpecMoments[i]);
        }
        res.put("SM02", SpecMoments[1]);
        res.put("SM03", SpecMoments[2]);
        res.put("SM04", SpecMoments[3]);
        res.put("SM05", SpecMoments[4]);
        res.put("SM06", SpecMoments[5]);
        res.put("SM07", SpecMoments[6]);
        res.put("SM08", SpecMoments[7]);


        return res;
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
}
