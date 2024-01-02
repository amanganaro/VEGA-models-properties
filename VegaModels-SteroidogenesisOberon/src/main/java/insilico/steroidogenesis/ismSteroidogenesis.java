package insilico.steroidogenesis;

import insilico.core.ad.item.ADIndexACF;
import insilico.core.ad.item.ADIndexAccuracy;
import insilico.core.ad.item.ADIndexConcordance;
import insilico.core.ad.item.ADIndexSimilarity;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.steroidogenesis.descriptors.EmbeddedDescriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.net.URL;

public class ismSteroidogenesis extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSteroidogenesis.class);
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_steroidogenesis.xml";

    private KNN knn = null;

    public ismSteroidogenesis() throws InitFailureException {
        super(ModelData);

        // Build model object
        URL KnnData = this.getClass().getResource("/data/knn_data.txt");
        knn = new KNN(KnnData);

        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted steroidogenesis activity";

        this.DescriptorsSize = 10;
        this.DescriptorsNames = new String[DescriptorsSize];


        DescriptorsNames[0] = "SRW9";
        DescriptorsNames[1] = "PubchemFP20";
        DescriptorsNames[2] = "PubchemFP37";
        DescriptorsNames[3] = "PubchemFP183";
        DescriptorsNames[4] = "PubchemFP189";
        DescriptorsNames[5] = "PubchemFP341";
        DescriptorsNames[6] = "PubchemFP342";
        DescriptorsNames[7] = "PubchemFP379";
        DescriptorsNames[8] = "PubchemFP418";
        DescriptorsNames[9] = "PubchemFP755";

        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

    }

    public ismSteroidogenesis(String ModelData) throws InitFailureException {
        super(ModelData);
    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.SRW9;
            Descriptors[1] = embeddedDescriptors.PubchemFP20;
            Descriptors[2] = embeddedDescriptors.PubchemFP37;
            Descriptors[3] = embeddedDescriptors.PubchemFP183;
            Descriptors[4] = embeddedDescriptors.PubchemFP189;
            Descriptors[5] = embeddedDescriptors.PubchemFP341;
            Descriptors[6] = embeddedDescriptors.PubchemFP342;
            Descriptors[7] = embeddedDescriptors.PubchemFP379;
            Descriptors[8] = embeddedDescriptors.PubchemFP418;
            Descriptors[9] = embeddedDescriptors.PubchemFP755;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        int prediction;
        prediction = knn.getPrediction(Descriptors);

        CurOutput.setMainResultValue(prediction);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(knn.getPrediction(Descriptors));
        } catch (Throwable ex) {
            log.warn("Unable to find label for Steroidogenesis value " + knn.getPrediction(Descriptors));
            Res[0] = Double.toString(knn.getPrediction(Descriptors));
        }
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;

    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
//        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
//        adq.setMoleculesForIndexSize(2);
//        if (!adq.Calculate(CurMolecule, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets threshold for AD indices
//        try {
//            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
//            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.6);
//            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.6);
//            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.6);
//        } catch (Throwable e) {
//            return InsilicoModel.AD_ERROR;
//        }
//
//        // Sets Range check
//        ADCheckDescriptorRange adrc = new ADCheckDescriptorRange();
//        if (!adrc.Calculate(TS, Descriptors, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets ACF check
//        ADCheckACF adacf = new ADCheckACF(TS);
//        if (!adacf.Calculate(CurMolecule, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets final AD index
//        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
//        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
//        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;
//
//        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
//        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
//                CurOutput.getADIndex(ADIndexConcordance.class),
//                CurOutput.getADIndex(ADIndexMaxError.class));
//        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
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
