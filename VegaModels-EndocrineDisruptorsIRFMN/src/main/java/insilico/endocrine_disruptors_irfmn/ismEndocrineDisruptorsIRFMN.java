package insilico.endocrine_disruptors_irfmn;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismEndocrineDisruptorsIRFMN extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismEndocrineDisruptorsIRFMN.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_endocrine_disruptors.xml";

    private final EndocrineECList ECList;
    private final EndocrineECSmarts ECSmarts;

    public ismEndocrineDisruptorsIRFMN() throws InitFailureException {
        super(ModelData);

        try {
            ECList = new EndocrineECList();
            ECSmarts = new EndocrineECSmarts();
        } catch (Exception e) {
            throw new InitFailureException("Unable to init EC list and smarts - " + e.getMessage());
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted ED activity";
        this.ResultsName[1] = "ED activity reason";

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

        int MainResult = 0;
        String MainResultText = "-";

        try {

            // EC List matching
            boolean OutECListFound = false;
            String OutECListAndSmarts = "";

                InsilicoMolecule ECListMolecule = ECList.Match(CurMolecule);
                if (ECListMolecule != null) {
                    MainResult = 1;
                    MainResultText = "Exact match with " + ECListMolecule.GetId();
                }

            // EC SMARTS matching
            if (MainResult != 1) {
                String res = ECSmarts.Match(CurMolecule);
                if (res != null) {
                    MainResult = 1;
                    MainResultText = "Match with category " + res;
                }
            }

        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(MainResult);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // ED classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for ED value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = MainResultText;

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;

    }

    @Override
    protected short CalculateAD() {

        return InsilicoModel.AD_ERROR;

    }

    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        if (CurOutput.getMainResultValue() == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (CurOutput.getMainResultValue() == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
