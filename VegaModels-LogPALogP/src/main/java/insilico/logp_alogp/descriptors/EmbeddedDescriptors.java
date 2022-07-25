package insilico.logp_alogp.descriptors;

import insilico.logp_alogp.descriptors.weights.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.interfaces.IAtomContainer;


@Log4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;

    private double[][] ConnAugMatrix;

    private double ALogP;

    public EmbeddedDescriptors(InsilicoMolecule Mol) {
        ALogP = MISSING_VALUE;
        CalculateDescriptors(Mol);
    }

    public void CalculateDescriptors(InsilicoMolecule Mol) {
        CalculateALogP(Mol);
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        this.setALogP(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setALogP(MISSING_VALUE);
            log.warn(e.getMessage());
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);

        double[] Frags = acf.GetAllValues();

        for (double d : Frags)
            if (d == MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            ALogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }
        this.setALogP(ALogP);
    }





}
