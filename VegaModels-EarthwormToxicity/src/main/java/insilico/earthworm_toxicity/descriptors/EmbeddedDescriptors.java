package insilico.earthworm_toxicity.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.molecule.tools.Manipulator;
import insilico.earthworm_toxicity.descriptors.weights.WeightsIState;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;


import java.util.ArrayList;

public class EmbeddedDescriptors {
	private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);
	private final int MISSING_VALUE = -999;
	public double B06_N_O = MISSING_VALUE;
	public double SaasC = MISSING_VALUE;
	public double B06_N_F = MISSING_VALUE;
	public double GATS1m = MISSING_VALUE;
	public double VE1signB_e = MISSING_VALUE;
	public double B03_N_Cl = MISSING_VALUE;
	public double B08_N_Cl = MISSING_VALUE;
	public double S_107 = MISSING_VALUE;

	public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
		InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION=true;
		InsilicoMolecule mol_norm = (InsilicoMolecule) mol.Clone();
		IAtomContainer ac = InsilicoMoleculeNormalization.Normalize(mol_norm.GetStructure());
		String SMI = SmilesMolecule.GenerateSmiles(ac);
		mol_norm.SetSMILESAndStructure(SMI, ac);
		InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION=false;
		CalculateAllDescriptors(mol_norm);
	}

	private void CalculateAllDescriptors(InsilicoMolecule mol) throws GenericFailureException, InvalidMoleculeException {
		CalculateBX(mol);
		CalculateSaasC(mol);
		CalculateGATS1m(mol);
		CalculateVE1signB_e(mol);
		CalculateS_107(mol);
	}


	private void CalculateVE1signB_e(InsilicoMolecule mol) throws InvalidMoleculeException, GenericFailureException {

		VE1signB_e = 0;

		int nSK = mol.GetStructure().getAtomCount();
		double[][] Mat = new double[nSK][nSK];
		try {
			Mat = mol.GetMatrixBurden();
			double[] w = new WeightsElectronegativity().getScaledWeights(mol.GetStructure());

			for (double val : w)
				if (val == MISSING_VALUE) {
					log.warn("unable_calculate_all_weights");
					VE1signB_e = MISSING_VALUE;
				}

			for (int i = 0; i < nSK; i++)
				Mat[i][i] = w[i];

		} catch (Exception e) {
			log.warn(e.getMessage());
			VE1signB_e = MISSING_VALUE;
		}

		Matrix DataMatrix = new Matrix(Mat);
		double[] eigenvalues;
		double[][] eigenVectors;

		try {
			EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
			eigenvalues = ed.getRealEigenvalues();
			eigenVectors = ed.getV().getArray();
		} catch (Throwable e) {
			log.warn("unable_eigenvalue" + e.getMessage());
			return;
		}

		int idxMin = 0;
		for (int idx=0; idx < eigenvalues.length; idx++){
			if (eigenvalues[idx] < eigenvalues[idxMin])
				idxMin = idx;
		}

		for (double[] eigenVector : eigenVectors) {
			VE1signB_e += eigenVector[idxMin];
		}
		VE1signB_e = Math.abs(VE1signB_e);
	}
	private void CalculateS_107(InsilicoMolecule mol) throws InvalidMoleculeException, GenericFailureException {
		S_107 = 0;

		int n = mol.GetStructure().getAtomCount();
		double[][] ConnMatrix = mol.GetMatrixConnectionAugmented();


		ArrayList VisitedS = new ArrayList();

		for(int i = 0; i < n ; i++) {
			if (mol.GetStructure().getAtom(i).getSymbol().equalsIgnoreCase("S")) {

				int ConnCount=0, R=0, S_SR=0;

				VisitedS.add(i);

				for (int ConnAtomIdx=0; ConnAtomIdx<n; ConnAtomIdx++){
					if (i==ConnAtomIdx) continue;
					if (ConnMatrix[i][ConnAtomIdx]>0) {
						ConnCount++;
						if (mol.GetStructure().getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("S")) {
							boolean AlreadyVisited = false;
							for (int w=0; w<VisitedS.size(); w++)
								if ((Integer) VisitedS.get(w) ==ConnAtomIdx) {
									AlreadyVisited = true;
									break;
								}
							if (!(AlreadyVisited)) {
								VisitedS.add(ConnAtomIdx);
								for (int k=0; k<n; k++){
									if (k==ConnAtomIdx) continue;
									if (ConnMatrix[k][ConnAtomIdx]>0) {
										S_SR++;
										break;
									}
								}
							}
						} else {
							R++;
						}

					}
				}

				if (ConnCount==2) {
					// R2S
					if (R==2) {
						S_107++;
						continue;
					}

					// RS-SR
					if (S_SR==1) {
						S_107++;
					}
				}
			}
		}
	}
	private void CalculateBX(InsilicoMolecule mol){

		String[][] ATOM_COUPLES = new String[3][2];
		ATOM_COUPLES[0] = new String[] {"N","O"};
		ATOM_COUPLES[1] = new String[] {"N","F"};
		ATOM_COUPLES[2] = new String[] {"N","Cl"};

		IAtomContainer m;
		try {
			m = mol.GetStructure();
		} catch (InvalidMoleculeException e) {
			B06_N_O = MISSING_VALUE;
			B03_N_Cl = MISSING_VALUE;
			B06_N_F = MISSING_VALUE;
			B08_N_Cl = MISSING_VALUE;
			return;
		}

		int nSK = m.getAtomCount();

		// Gets matrices
		int[][] TopoMat = null;
		try {
			TopoMat = mol.GetMatrixTopologicalDistance();
		} catch (GenericFailureException e) {
			log.warn("Invalid structure, unable to calculate: " + mol);
			B06_N_O = MISSING_VALUE;
			B03_N_Cl = MISSING_VALUE;
			B06_N_F = MISSING_VALUE;
			B08_N_Cl = MISSING_VALUE;
			return;
		}

		B06_N_O = 0;
		B03_N_Cl = 0;
		B06_N_F = 0;
		B08_N_Cl = 0;

		for (int d = 0; d < ATOM_COUPLES.length; d++) {
			for (int i=0; i<nSK; i++) {
				if (m.getAtom(i).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][0])) {
					for (int j=0; j<nSK; j++) {
						if (i==j) continue;
						if (m.getAtom(j).getSymbol().equalsIgnoreCase(ATOM_COUPLES[d][1])) {
							if ((d == 0) && TopoMat[i][j]==6)
								B06_N_O = 1;
							if ((d == 1) && TopoMat[i][j]==6)
								B06_N_F = 1;
							if ((d == 2) && TopoMat[i][j]==3)
								B03_N_Cl = 1;
							if ((d == 2) && TopoMat[i][j]==8)
								B08_N_Cl = 1;
						}

					}
				}
			}
		}
	}
	private void CalculateSaasC(InsilicoMolecule mol){
		SaasC = 0;

		IAtomContainer m;
		try {
			m = mol.GetStructure();
		} catch (InvalidMoleculeException e) {
			log.warn("invalid_structure");
			SaasC = MISSING_VALUE;
			return;
		}

		int[][] TopoMatrix;
		try {
			TopoMatrix =  mol.GetMatrixTopologicalDistance();
		} catch (Exception e) {
			log.warn(e.getMessage());
			SaasC = MISSING_VALUE;
			return;
		}

		int nSK = m.getAtomCount();

		// Get I-States weights
		WeightsIState w_istate = new WeightsIState();
		double[] w_is = w_istate.getWeights(m, false);
		for (double val : w_is)
			if (val == MISSING_VALUE) {
				log.warn("unable_calculate_all_estates");
				SaasC = MISSING_VALUE;
				return;
			}

		// Calculate E-States
		double[] w_es = new double[nSK];
		for (int at = 0; at<nSK; at++) {
			double sumDeltaI = 0;

			for (int j = 0; j < nSK; j++)
				if (at != j)
					sumDeltaI += (w_is[at] - w_is[j]) / Math.pow((double) TopoMatrix[at][j] + 1.0, 2.0);

			w_es[at] = w_is[at] + sumDeltaI;
		}

		for (int at=0; at<m.getAtomCount(); at++) {

			IAtom curAt = m.getAtom(at);
			if (!curAt.getSymbol().equals("C"))
				continue;

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

			// Groups count
			if ((nSng == 1) && (nDbl == 0) &&
					(nTri == 0) && (nAr == 2) &&
					(nH == 0) && (Charge == 0)) {
				SaasC += w_es[at];
			}
		}
	}
	private void CalculateGATS1m(InsilicoMolecule mol) {

			GATS1m = MISSING_VALUE;

			IAtomContainer m;
			try {
				IAtomContainer orig_m = mol.GetStructure();
				m = Manipulator.AddHydrogens(orig_m);
			} catch (InvalidMoleculeException | GenericFailureException e) {
				log.warn("invalid_structure");
				GATS1m = MISSING_VALUE;
				return;
			}

			// Gets matrices
			int[][] TopoMatrix;
			try {
				TopoMatrix = TopoDistanceMatrix.getMatrix(m);
			} catch (Exception e) {
				log.warn(e.getMessage());
				GATS1m = MISSING_VALUE;
				return;
			}

			int nSK = m.getAtomCount();
			double[] w =  new WeightsMass().getScaledWeights(m);

			// If one or more weights are not available, sets all to missing value
			for (int i=0; i<nSK; i++)
				if (w[i] == MISSING_VALUE) {
					return;
				}

			// Calculates weights averages
			double wA = 0;
			for (int i=0; i<nSK; i++)
				wA += w[i];
			wA = wA / ((double) nSK);

			double GearyAC=0;
			double denom = 0, delta = 0;

			for (int i=0; i<nSK; i++) {

				denom += Math.pow((w[i] - wA), 2);

				for (int j=0; j<nSK; j++)
					if (TopoMatrix[i][j] == 1) {
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
			GATS1m = GearyAC;
	}
}

