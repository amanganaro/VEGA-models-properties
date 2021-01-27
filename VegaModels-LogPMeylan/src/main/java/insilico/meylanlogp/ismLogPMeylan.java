package insilico.meylanlogp;

import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.meylanlogp.descriptors.ADIndexADILogP;
import insilico.meylanlogp.descriptors.EmbeddedDescriptors;


public class ismLogPMeylan extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_logp_meylan.xml";

    private boolean SkipExperimentalInAD;

    public ismLogPMeylan()
            throws InitFailureException {
        super(ModelData);

        SkipExperimentalInAD = false;

        // Define no. of descriptors
        this.DescriptorsSize = 1;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "MeylanLogP";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LogP";

        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();

    }

    public void SetSkipExperimentalInAD(boolean skip) {
        this.SkipExperimentalInAD = skip;
    }

//    @Override
//    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
//
//        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
//        DescriptorBlock desc;
//
//        // Meylan LogP
//        desc = new MeylanLogP();
//        blocks.add(desc);
//
//        return blocks;
//    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors descriptors = new EmbeddedDescriptors();
            descriptors.CalculateAllDescriptors(CurMolecule);


            Descriptors[0] = descriptors.LogP;
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }



    @Override
    protected short CalculateModel() {

        CurOutput.setMainResultValue(Descriptors[0]);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Descriptors[0])); // LogP
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setSkipExperimental(SkipExperimentalInAD);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.9, 0.75);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.5);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets final AD index
        ADIndexADI ADI = new ADIndexADILogP(adq.getIndexADI(), CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));

        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

        double LogP_threshold_green = 3;
        double LogP_threshold_red = 8;

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val < LogP_threshold_green)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val < LogP_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
    }

}
