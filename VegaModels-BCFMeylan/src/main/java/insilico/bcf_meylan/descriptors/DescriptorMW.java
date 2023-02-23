package insilico.bcf_meylan.descriptors;

import insilico.bcf_meylan.descriptors.weights.Mass;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;


public class DescriptorMW {

    public double Mw = -999;
    private static final Logger log = LogManager.getLogger(DescriptorMW.class);


    public DescriptorMW(InsilicoMolecule Mol){
        CalculateMw(Mol);
    }

    private void CalculateMw(InsilicoMolecule Mol){

        Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            Mw = -999;
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[] H = new int[nSK];

            for (int i=0; i<nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) {
                    log.warn(e.getMessage());
                }
            }

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");

            for (int i=0; i<nSK; i++) {
                if (wMass[i] == -999)
                    Mw = -999;
            }

            for (int i=0; i<nSK; i++) {
                if (Mw != -999) {
                    Mw += wMass[i];
                    if (H[i] > 0) {
                        Mw += HMass * H[i];
                    }
                }
            }

        } catch (Throwable e) {
            log.warn(e.getMessage());
            Mw = -999;
        }
    }

}
