package insilico.ppb_coral;

import insilico.core.ad.item.ADIndexACF;
import insilico.core.ad.item.ADIndexAccuracy;
import insilico.core.ad.item.ADIndexConcordance;
import insilico.core.ad.item.ADIndexSimilarity;
import insilico.core.coral.CoralModel;
import insilico.core.coral.models.PPB.CoralPPB;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;

@Slf4j
public class ismPPBCoral extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_ppb_coral.xml";

    private CoralModel CoralPPB;

    public ismPPBCoral() throws InitFailureException {
        super(ModelData);
        try {
            CoralPPB = new CoralPPB();
        } catch (Exception ex){
            throw new InitFailureException("Unable to init coral model");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "PPB [square root of fraction unbound]";

        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

    }

    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
        return new ArrayList<>();
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
        double Prediction;
        try {
            Prediction = CoralPPB.Predict(this.CurMolecule.getInputSMILES());
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        this.ResultsSize = 1;
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Prediction);


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

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == -1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
