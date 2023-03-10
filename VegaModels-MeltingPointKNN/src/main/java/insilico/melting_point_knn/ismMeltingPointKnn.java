package insilico.melting_point_knn;

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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ismMeltingPointKnn extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismMeltingPointKnn.class);

    private static final long serialVersionUID = 1L;
    private insilicoKnnQuantitative KNN1;
//    private insilicoKnnQuantitative KNN3;
    private insilicoKnnPrediction Knn1Prediction;
//    private insilicoKnnPrediction Knn3Prediction;

    private static final String ModelData = "/data/model_melting_point_knn.xml";

    public ismMeltingPointKnn()
                throws InitFailureException {
            super(ModelData);

            // Build model object
            KNN1 =  new insilicoKnnQuantitative();
            KNN1.setNeighboursNumber(4);
            KNN1.setMinSimilarity(0.7);
            KNN1.setMinSimilarityForSingleResult(0.75);
            KNN1.setEnhanceWeightFactor(3);
            KNN1.setUseExperimentalRange(false);

           // Build model object
//            KNN3 =  new insilicoKnnQuantitative();
//            KNN3.setNeighboursNumber(4);
//            KNN3.setMinSimilarity(0.8);
//            KNN3.setMinSimilarityForSingleResult(0.9);
//            KNN3.setEnhanceWeightFactor(3);
//            KNN3.setUseExperimentalRange(false);


           // Define no. of descriptors
            this.DescriptorsSize = 0;
            this.DescriptorsNames = new String[DescriptorsSize];

            // Defines results
            this.ResultsSize = 2;
            this.ResultsName = new String[ResultsSize];
            this.ResultsName[0] = "Predicted melting point [°C]";
            this.ResultsName[1] = "Molecules used for prediction";

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

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }


    @Override
    protected short CalculateModel() {
        try {
            Knn1Prediction = KNN1.Calculate(CurMolecule, TS, this.KnnSkipExperimental);
//            Knn3Prediction = KNN3.Calculate(CurMolecule, TS, this.KnnSkipExperimental);
//            BufferedWriter writer = new BufferedWriter(new FileWriter(new File("VegaModels-MeltingPointKNN/src/main/resources/smilesMP_out.txt"),true));
//            writer.write(CurMolecule.GetSMILES() + "\t" + Knn1Prediction.getPrediction() + "\t" + Knn3Prediction.getPrediction() + "\t");
//            writer.flush();
//            System.out.print(CurMolecule.GetSMILES() + "\t" + Knn1Prediction.getPrediction() + "\t" + Knn3Prediction.getPrediction());
        } catch (GenericFailureException ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Knn1Prediction.getPrediction());

        String[] Res = new String[ResultsSize];
        if (Knn1Prediction.getPrediction() == Descriptor.MISSING_VALUE)
            Res[0] = "-";
        else
            Res[0] = String.valueOf(Format_2D.format(Knn1Prediction.getPrediction()));
        Res[1] = String.valueOf(Knn1Prediction.getNeighbours().size());
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        // Ad is performed on the K molecules used for the KNN model
        adq.setMoleculesForIndexSize(Knn1Prediction.getNeighbours().size());

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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(50.0, 10.0);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(50.0, 10.0);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(50.0, 10.0);
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

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        if (CurOutput.getMainResultValue() == Descriptor.MISSING_VALUE) {

            CurOutput.setAssessment("N/A");
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_NA, "N/A"));

        } else {

            ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0], "°C");

        }

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this);
        TSK.SerializeToFile(DatName);
    }

}
