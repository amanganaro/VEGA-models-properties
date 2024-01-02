package insilico.steroidogenesis.descriptors;

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
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.earthworm_toxicity.descriptors.weights.WeightsIState;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;


import java.util.ArrayList;

public class EmbeddedDescriptors {
	private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);
	private final int MISSING_VALUE = -999;

	public double SRW9 = MISSING_VALUE;
	public double PubchemFP20 = MISSING_VALUE;
	public double PubchemFP37 = MISSING_VALUE;
	public double PubchemFP183 = MISSING_VALUE;
	public double PubchemFP189 = MISSING_VALUE;
	public double PubchemFP341 = MISSING_VALUE;
	public double PubchemFP342 = MISSING_VALUE;
	public double PubchemFP379 = MISSING_VALUE;
	public double PubchemFP418 = MISSING_VALUE;
	public double PubchemFP755 = MISSING_VALUE;

	public double prova = MISSING_VALUE;

	public EmbeddedDescriptors(InsilicoMolecule mol) throws Exception {
		CalculateAllDescriptors(mol);
	}

	private void CalculateAllDescriptors(InsilicoMolecule mol) throws GenericFailureException, InvalidMoleculeException, CDKException {
		CalculateSRW9(mol);

		// fixed - SRW was normalized in the original model from Marco, should be plain range scaling on max val = 8.65364531455174
		SRW9 /= 8.65364531455174;

		CalculatePubChemFPs(mol);
	}

	private void CalculatePubChemFPs(InsilicoMolecule mol) throws GenericFailureException, CDKException {
		IAtomContainer m;

		try {
			m = mol.GetStructure();
		} catch (InvalidMoleculeException e) {
			log.warn("Invalid structure, unable to calculate paths");
			throw new GenericFailureException("invalid structure");
		}

		PubchemFingerprinterVega.CountRings countRings = new PubchemFingerprinterVega.CountRings(m);
		PubchemFingerprinterVega.CountElements countElements = new PubchemFingerprinterVega.CountElements(m);
		PubchemFingerprinterVega.CountSubstructures countSubstructures = new PubchemFingerprinterVega.CountSubstructures(m);

		PubchemFP20 = (countElements.getCount("O") >= 4) ? 1 : 0;
		PubchemFP37 = (countElements.getCount("Cl") >= 1) ? 1 : 0;

		// fixed - this FP is wrong but Padel calculates it like this
		int nRingNUnsat = countRings.countUnsaturatedNitrogenContainingRing(6);
		int nRingNArom =  countRings.countAromaticNitrogenContainingRing(6);
		PubchemFP183 = (nRingNUnsat+nRingNArom >= 1) ? 1 : 0;
//		PubchemFP183 = (countRings.countUnsaturatedNitrogenContainingRing(6) >= 1) ? 1 : 0;

		// fixed - this FP is wrong but Padel calculates it like this
		int nOnlyCarbUnsat = countRings.countUnsaturatedCarbonOnlyRing(6);
		int nBenzene = countSubstructures.countSubstructure("c1ccccc1");
		PubchemFP189 = ( nOnlyCarbUnsat+nBenzene >= 2) ? 1 : 0;
//		PubchemFP189 = (countRings.countUnsaturatedCarbonOnlyRing(6) >= 2) ? 1 : 0;

		PubchemFP341 = (countSubstructures.countSubstructure("[#6](~[#6])(~[#6])(~[#8])") > 0) ? 1 : 0;
		PubchemFP342 = (countSubstructures.countSubstructure("[#6](~[#6])(~[Cl])") > 0) ? 1 : 0;
		PubchemFP379 = (countSubstructures.countSubstructure("[#6](~[#7])(:n)") > 0) ? 1 : 0;
		PubchemFP418 = (countSubstructures.countSubstructure("[#6]=,:[#7]") > 0) ? 1 : 0;
		PubchemFP755 = (countSubstructures.countSubstructure("[#6]c1c([#6])cccc1") > 0) ? 1 : 0;

		prova = countSubstructures.countSubstructure("c1ccccc1");

	}
	private void CalculateSRW9(InsilicoMolecule mol) throws GenericFailureException {

		IAtomContainer m;

		try {
			m = mol.GetStructure();
		} catch (InvalidMoleculeException e) {
			log.warn("Invalid structure, unable to calculate paths");
			throw new GenericFailureException("invalid structure");
		}

		int[][] AdjMat;
		try {
			AdjMat = mol.GetMatrixAdjacency();
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

		int nSK = m.getAtomCount();

		double[] Self_Returning_Walk_Counts = new double[10];

		double[][] AdjPow = new double[nSK][nSK];

		Self_Returning_Walk_Counts[0] = nSK;

		int Pow_k = 1;
		while (Pow_k < 9) {

			for (int i=0; i<nSK; i++) {
				for (int j = 0; j < nSK; j++) {
					double ww = 0;
					for (int k = 0; k < nSK; k++)
						ww += AdjMatDbl[i][k] * AdjxAdj[k][j];
					AdjPow[i][j] = ww;
				}
			}

			for (int i=0; i<nSK; i++)
				Self_Returning_Walk_Counts[Pow_k] = Self_Returning_Walk_Counts[Pow_k] + AdjPow[i][i];

			Pow_k++;

			for (int i=0; i<nSK; i++)
				for (int j = 0; j < nSK; j++)
					AdjxAdj[i][j] = AdjPow[i][j];

		}

		SRW9 = Math.log(1 + Self_Returning_Walk_Counts[8]);

	}

}

