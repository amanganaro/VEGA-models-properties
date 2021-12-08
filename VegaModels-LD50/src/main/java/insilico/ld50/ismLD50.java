package insilico.ld50;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.knn.insilicoKnnQuantitative;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;

public class ismLD50 extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_ld50_knn.xml";

    private final insilicoKnnQuantitative KNN;
    private insilicoKnnPrediction KnnPrediction;
    private double MW;


    public ismLD50()
            throws InitFailureException {
        super(ModelData);

        // Build model object
        KNN =  new insilicoKnnQuantitative();
        KNN.setNeighboursNumber(3);
        KNN.setMinSimilarity(0.8);
        KNN.setMinSimilarityForSingleResult(0.85);
        KNN.setEnhanceWeightFactor(3);
        KNN.setExperimentalRange(2);
        KNN.setUseExperimentalRange(true);

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted log LD50 [log(mmol/Kg)]";
        this.ResultsName[1] = "Predicted log LD50 [mg/Kg]";
        this.ResultsName[2] = "Molecules used for prediction";
        this.ResultsName[3] = "Experimental value [mg/Kg]";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();

    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            // MW in constitutional is given as a SCALED
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            MW = CarbonWeight * CurMolecule.GetBasicDescriptorByName("MW").getValue();

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        try {
            KnnPrediction = KNN.Calculate(CurMolecule, TS);
        } catch (GenericFailureException ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(KnnPrediction.getPrediction());

        String[] Res = new String[ResultsSize];
        if (KnnPrediction.getPrediction() == Descriptor.MISSING_VALUE) {
            Res[0] = "-";
            Res[1] = "-";
        } else {
            Res[0] = String.valueOf(Format_3D.format(KnnPrediction.getPrediction())); // log(mmol/kg)
            double ConvertedValue = Math.pow(10, KnnPrediction.getPrediction()) * MW;
            if (ConvertedValue>1)
                Res[1] = Format_2D.format(ConvertedValue); // mg/kg
            else
                Res[1] = Format_4D.format(ConvertedValue); // mg/kg
        }
        Res[2] = String.valueOf(KnnPrediction.getNeighbours().size());
        Res[3] = "-";
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {
        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        // Ad is performed on the K molecules used for the KNN model
        adq.setMoleculesForIndexSize(KnnPrediction.getNeighbours().size());

        // (only retrieve similar molecules if n.a. prediction)
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == Descriptor.MISSING_VALUE) {
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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.75, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.2, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.2, 0.6);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.2, 0.6);
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

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.75, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);

        // Add transformed (mg/L) experimental if needed
        if (CurOutput.HasExperimental()) {
            double ConvertedValue = Math.pow(10, CurOutput.getExperimental()) * MW;
            if (ConvertedValue>1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedValue); // mg/kg
            else
                CurOutput.getResults()[3] = Format_4D.format(ConvertedValue); // mg/kg
        }

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {

        double LC_threshold_red = 1; // in mg/kg
        double LC_threshold_orange = 10; // in mg/kg
        double LC_threshold_yellow = 100; // in mg/kg

        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set (mg/kg) if available

        if (CurOutput.getMainResultValue() == Descriptor.MISSING_VALUE) {

            CurOutput.setAssessment("N/A");
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_NA, "N/A"));
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

        } else {

            String ADItemWarnings =
                    ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());

            String Result = CurOutput.getResults()[1] + " mg/kg";

            switch (CurOutput.getADI().GetAssessmentClass()) {
                case ADIndex.INDEX_LOW:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_LOW, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_LOW, Result, ADItemWarnings));
                    break;
                case ADIndex.INDEX_MEDIUM:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_MEDIUM, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_MEDIUM, Result, ADItemWarnings));
                    break;
                case ADIndex.INDEX_HIGH:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_HIGH, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_HIGH, Result));
                    if (!ADItemWarnings.isEmpty())
                        CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() +
                                String.format(MessagesAD.ASSESS_LONG_ADD_ISSUES, ADItemWarnings));
                    break;
            }

            // Override assessment if experimental value is available
            if (CurOutput.HasExperimental()) {
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/Kg", CurOutput.getAssessment()));
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/Kg"));
            }


            // Sets assessment status
            double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
            Val = Math.pow(10, Val) * MW;
            if (Val < LC_threshold_red)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
            else if (Val < LC_threshold_orange)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
            else if (Val < LC_threshold_yellow)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
            else
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        }
    }


    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this, true, false);
        TSK.SerializeToFile(DatName);
    }

}
