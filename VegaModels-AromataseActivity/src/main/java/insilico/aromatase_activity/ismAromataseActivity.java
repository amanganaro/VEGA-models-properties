package insilico.aromatase_activity;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

@Slf4j
public class ismAromataseActivity extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_aromatase_tox21.xml";

    private AromataseActivitySMARTS SAs;

    public ismAromataseActivity() throws InitFailureException {
        super(ModelData);

        try {
            SAs = new AromataseActivitySMARTS();
        } catch (Exception e) {
            throw new InitFailureException("Unable to init smarts - " + e.getMessage());
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted activity";
        this.ResultsName[1] = "Alerts for activity";
        this.ResultsName[2] = "Alerts for inactivity";
        this.ResultsName[3] = "Alerts for agonist activity";
        this.ResultsName[4] = "Alerts for antagonist activity";

        // Define AD items
        this.ADItemsName = new String[0];
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        int MainResult;

        try {

            SAs.Match(CurMolecule);

        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        if (!SAs.Matches_Active.isEmpty()) {
            // Active
            if (!SAs.Matches_Active_Agonist.isEmpty()) {
                MainResult = 1;
            } else if (!SAs.Matches_Active_Antagonist.isEmpty()) {
                MainResult = 2;
            } else MainResult = 3;
        } else if (!SAs.Matches_Inactive.isEmpty()) {
            // Inactive
            MainResult = 0;
        } else {
            // not predicted
            MainResult = -1;
        }

        CurOutput.setMainResultValue(MainResult);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // aromatase classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for ED value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = AromataseActivitySMARTS.FormatAlertArray(SAs.Matches_Active);
        Res[2] = AromataseActivitySMARTS.FormatAlertArray(SAs.Matches_Inactive);
        Res[3] = AromataseActivitySMARTS.FormatAlertArray(SAs.Matches_Active_Agonist);
        Res[4] = AromataseActivitySMARTS.FormatAlertArray(SAs.Matches_Active_Antagonist);

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        if (CurOutput.getMainResultValue() == -1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else if (CurOutput.getMainResultValue() == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);

    }
}
