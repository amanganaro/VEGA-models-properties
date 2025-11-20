package insilico.ttr;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import insilico.ttr.model.RandomForestPredictor;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismTTR extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismTTR.class);
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_ttr.xml";

    private final RandomForestPredictor model_rf;


    public ismTTR() throws InitFailureException {
        super(ModelData);

        // Build model object
        model_rf = new RandomForestPredictor();

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted TTR activity";
        this.ResultsName[1] = "RF model confidence";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

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

        try {

            RandomForestPredictor.Prediction pred = model_rf.predict(CurMolecule.getInputSMILES());
            int final_res = -1;
            try {
                final_res = Integer.parseInt(pred.prediction);
            } catch (NumberFormatException e) {
                // do nothing, prediction alredy set to -1
            }

            CurOutput.setMainResultValue(final_res);

            String[] Res = new String[ResultsSize];
            try {
                Res[0] = this.GetTrainingSet().getClassLabel(final_res);
            } catch (Throwable ex) {
                log.warn("Unable to find label for TTR value " + final_res);
                Res[0] = Double.toString(final_res);
            }
            Res[1] = "" + pred.confidence;

            CurOutput.setResults(Res);

        } catch (Exception e) {
            return MODEL_ERROR;
        }

        return MODEL_CALCULATED;

    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(0);
        adq.AddMappingToNegativeValue(1);
        adq.setMoleculesForIndexSize(2);

        // (only retrieve similar molecules if n.a. prediction)
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == -1) {
            try {
                adq.SetSimilarMolecules(CurMolecule, CurOutput);
            } catch (GenericFailureException ex) {
                // do nothing
            }
            return InsilicoModel.AD_ERROR;
        }

        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.85, 0.7);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.85, 0.7);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.85, 0.7);
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(true);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this);
        TSK.SerializeToFile(DatName);
    }
}
