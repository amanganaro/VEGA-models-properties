package insilico.meylanlogp.descriptors;

import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAMeylanLogPAdditionalFragments;
import insilico.core.alerts.builders.SAMeylanLogPCorrectionFragments;
import insilico.core.alerts.builders.SAMeylanLogPFragments;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.interfaces.IAtomContainer;


public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private final int MISSING_VALUE = -999;
    public double LogP;


    public EmbeddedDescriptors(){
        LogP = MISSING_VALUE;
    }

    public void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateMeylanLogP(Mol);
    }


    private void CalculateMeylanLogP(InsilicoMolecule mol) {
        LogP = 0;

        try {
            IAtomContainer m = mol.GetStructure();

            SAMeylanLogPFragments Frags = new SAMeylanLogPFragments();
            AlertList f_base = Frags.Calculate(mol);
            SAMeylanLogPAdditionalFragments FragsAdd = new SAMeylanLogPAdditionalFragments();
            AlertList f_add = FragsAdd.Calculate(mol);
            SAMeylanLogPCorrectionFragments FragsCorr = new SAMeylanLogPCorrectionFragments();
            AlertList f_corr = FragsCorr.Calculate(mol);

            double Coefficient = Frags.GetCoefficient() + FragsAdd.GetCoefficient();
            double Correction = FragsCorr.GetCoefficient();

            // main eq
            LogP = 0.2290 + Coefficient + Correction;

            // lower bound
            if (LogP < -5.0)
                LogP = -5.0;

        } catch (InitFailureException | InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            LogP = MISSING_VALUE;
        }
    }


}
