package insilico.steroidogenesisedscreen;

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
import insilico.steroidogenesisedscreen.model.SteroidogenesisConsModel;
import insilico.steroidogenesisedscreen.model.SteroidogenesisConsResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismSteroidogenesisEDScreen extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSteroidogenesisEDScreen.class);
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_steroidogenesis_eds.xml";

    private final SteroidogenesisConsModel consModel;


    public ismSteroidogenesisEDScreen() throws InitFailureException {
        super(ModelData);

        // Build model object
        consModel = new SteroidogenesisConsModel();

        this.ResultsSize = 7;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted steroidogenesis activity";
        this.ResultsName[1] = "KNN model predicted activity";
        this.ResultsName[2] = "KNN model - similarity";
        this.ResultsName[3] = "KNN model - n. of positives";
        this.ResultsName[4] = "KNN model - n. of negatives";
        this.ResultsName[5] = "RF model predicted activity";
        this.ResultsName[6] = "RF model - confidence";

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

        SteroidogenesisConsResult res = consModel.predict(CurMolecule);

        CurOutput.setMainResultValue(res.Result);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(res.Result);
        } catch (Throwable ex) {
            log.warn("Unable to find label for Steroidogenesis value " + res.Result);
            Res[0] = Double.toString(res.Result);
        }
        try {
            Res[1] = this.GetTrainingSet().getClassLabel(res.KNN_result);
        } catch (Throwable ex) {
            log.warn("Unable to find label for Steroidogenesis value " + res.KNN_result);
            Res[1] = Double.toString(res.KNN_result);
        }
        Res[2] = "" + res.KNN_similarity;
        Res[3] = "" + res.KNN_n_positive;
        Res[4] = "" + res.KNN_n_negative;
        try {
            Res[5] = this.GetTrainingSet().getClassLabel(res.RF_result);
        } catch (Throwable ex) {
            log.warn("Unable to find label for Steroidogenesis value " + res.RF_result);
            Res[5] = Double.toString(res.RF_result);
        }
        Res[6] = "" + res.RF_confidence;
        CurOutput.setResults(Res);

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
