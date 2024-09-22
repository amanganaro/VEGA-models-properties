package insilico.steroidogenesis;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismSteroidogenesis extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSteroidogenesis.class);
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_steroidogenesis.xml";

//    private KNN knn = null;
    private final FpKNN FingerprintKNN;


    public ismSteroidogenesis() throws InitFailureException {
        super(ModelData);

        // Build model object
//        URL KnnData = this.getClass().getResource("/data/knn_data_only_train.txt");
//        knn = new KNN(KnnData);
        FingerprintKNN = new FpKNN();

        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted steroidogenesis activity";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

//        this.DescriptorsSize = 10;
//        this.DescriptorsNames = new String[DescriptorsSize];
//
//        DescriptorsNames[0] = "SRW9";
//        DescriptorsNames[1] = "PubchemFP20";
//        DescriptorsNames[2] = "PubchemFP37";
//        DescriptorsNames[3] = "PubchemFP183";
//        DescriptorsNames[4] = "PubchemFP189";
//        DescriptorsNames[5] = "PubchemFP341";
//        DescriptorsNames[6] = "PubchemFP342";
//        DescriptorsNames[7] = "PubchemFP379";
//        DescriptorsNames[8] = "PubchemFP418";
//        DescriptorsNames[9] = "PubchemFP755";

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

//            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
//
//            Descriptors[0] = embeddedDescriptors.SRW9;
//            Descriptors[1] = embeddedDescriptors.PubchemFP20;
//            Descriptors[2] = embeddedDescriptors.PubchemFP37;
//            Descriptors[3] = embeddedDescriptors.PubchemFP183;
//            Descriptors[4] = embeddedDescriptors.PubchemFP189;
//            Descriptors[5] = embeddedDescriptors.PubchemFP341;
//            Descriptors[6] = embeddedDescriptors.PubchemFP342;
//            Descriptors[7] = embeddedDescriptors.PubchemFP379;
//            Descriptors[8] = embeddedDescriptors.PubchemFP418;
//            Descriptors[9] = embeddedDescriptors.PubchemFP755;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        double prediction = -1;
        try {
            prediction = FingerprintKNN.Calculate(CurMolecule, this.KnnSkipExperimental);
        } catch (Exception e) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(prediction);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(prediction);
        } catch (Throwable ex) {
            log.warn("Unable to find label for Steroidogenesis value " + prediction);
            Res[0] = Double.toString(prediction);
        }
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

//    @Override
//    public void ProcessTrainingSet() throws Exception {
//        this.setSkipADandTSLoading(false);
//        TrainingSet TSK = new TrainingSet();
//        String TSPath = this.getInfo().getTrainingSetURL();
//        String[] buf = TSPath.split("/");
//        String DatName = buf[buf.length-1];
//        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
//        TSK.Build(TSPath, this);
//        TSK.SerializeToFile(DatName);
//    }
}
