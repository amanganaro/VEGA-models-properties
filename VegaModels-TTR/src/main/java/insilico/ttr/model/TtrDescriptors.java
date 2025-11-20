package insilico.ttr.model;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

public class TtrDescriptors {

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

        int[][] TopoMatrix;
        try {
            TopoMatrix =  mol.GetMatrixTopologicalDistance();
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate topological distance matrix");
        }


        //// SsOH

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
        }

        double D02_SsOH = 0;
        for (int at=0; at<m.getAtomCount(); at++) {
            IAtom curAt = m.getAtom(at);
            if (curAt.getAtomicNumber() != 8)
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
            if ( (nBnd==1) && (nSng==1) && (nAr == 0) && (nH == 1) && (Charge == 0) )
                D02_SsOH += w_es[at];

        }



        //////// descriptors with H-filled structure

        IAtomContainer curMolH = Manipulator.AddHydrogens(m);

        int[][] TopoDistMatH;
        try {
            TopoDistMatH = TopoDistanceMatrix.getMatrix(curMolH);
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate matrix - " + e);
        }
        int nSKH = curMolH.getAtomCount();


        //// MATS1e

        double[] w = (new WeightsElectronegativity()).getScaledWeights(curMolH);

        for (int i=0; i<nSKH; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                throw new GenericFailureException("missing electronegativity weight");

        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSKH; i++)
            wA += w[i];
        wA = wA / ((double) nSKH);

        int lag = 1;

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSKH; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSKH; j++)
                if (TopoDistMatH[i][j] == lag) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0)
                MoranAC = 1;
            else
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSKH)) * denom);
        }

        double D01_MATS1e = MoranAC;



        double[] res = new double[2];
        res[0] = D01_MATS1e;
        res[1] = D02_SsOH;

        return res;
    }

    public final static String[] DESCRIPTOR_NAMES = {
            "MATS1e",
            "SsOH"
    };

    // values calculated from the saved TRAINING SET only
    private final double[] DESC_MAX =
            {0.310523, 36.330693};
    private final double[] DESC_MIN =
            {-0.301202, 0};

    public double[] Scale(double[] descriptors) {
        double[] scaled = new double[descriptors.length];
        for (int i=0; i<descriptors.length; i++)
            scaled[i] = ( descriptors[i] - DESC_MIN[i] ) / (DESC_MAX[i] - DESC_MIN[i]);
        return scaled;
    }

}
