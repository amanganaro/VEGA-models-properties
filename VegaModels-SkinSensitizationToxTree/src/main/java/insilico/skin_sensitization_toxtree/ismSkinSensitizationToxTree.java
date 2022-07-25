package insilico.skin_sensitization_toxtree;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;

@Log4j
public class ismSkinSensitizationToxTree extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization.xml";

    private final static boolean CALCULATE_AD = false;

    private ToxTreeSkinClassification ToxTreeSkin;


    public ismSkinSensitizationToxTree() throws InitFailureException {
        super(ModelData);

        ToxTreeSkin = new ToxTreeSkinClassification();

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted class";
        this.ResultsName[1] = "Predicted class description";

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
        ToxTreeSkinClassification.SkinClass oToxTree;

        try {
            oToxTree = ToxTreeSkin.Calculate(CurMolecule);
        } catch (GenericFailureException e) {
            return MODEL_ERROR;
        }


        CurOutput.setMainResultValue(oToxTree.getId());
        String[] Res = new String[ResultsSize];
        Res[0] = oToxTree.getName();
        Res[1] = oToxTree.getDescription();

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;


    }

    @Override
    protected short CalculateAD() {
        // for now - AD not calculated
        return InsilicoModel.AD_ERROR;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
