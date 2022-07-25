package insilico.zebrafish_coral.descriptor;

import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

@Log4j
@Data
public class EmbeddedDescriptor {

    private double MISSING_VALUE = -999;
    private double MW = MISSING_VALUE;

    public EmbeddedDescriptor(InsilicoMolecule Mol){
        CalculateMW(Mol);
    }

    private void CalculateMW(InsilicoMolecule Mol){
        this.setMW(0);
        double Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            this.setMW(MISSING_VALUE);
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
                    this.setMW(MISSING_VALUE);
            }

            for (int i=0; i<nSK; i++) {
                if (Mw != -999) {
                    Mw += wMass[i];
                    if (H[i] > 0) {
                        Mw += HMass * H[i];
                    }
                }
            }
            this.setMW(Mw);

        } catch (Throwable e) {
            log.warn(e.getMessage());
            this.setMW(MISSING_VALUE);
        }

    }


}
