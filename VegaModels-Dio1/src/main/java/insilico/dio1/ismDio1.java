package insilico.dio1;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import insilico.dio1.model.modelGBTpmml.GBTMolecularClassifier;
import insilico.dio1.model.modelGBTpmml.GBTResult;
import insilico.dio1.model.modelKNN.DioMolecularKNNClassifierFixed;
import insilico.dio1.model.modelSarpy.DioSarpy;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismDio1 extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismDio1.class);
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_dio1.xml";

    private final GBTMolecularClassifier model_gbt;
    private final DioMolecularKNNClassifierFixed model_knn;
    private final DioSarpy model_sarpy;


    public ismDio1() throws InitFailureException {
        super(ModelData);

        // Build model object
        model_gbt = new GBTMolecularClassifier();
        model_knn = new DioMolecularKNNClassifierFixed();
        model_sarpy = new DioSarpy();

        this.ResultsSize = 8;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted DIO1 activity";
        this.ResultsName[1] = "KNN model predicted activity";
        this.ResultsName[2] = "KNN model - similarity";
        this.ResultsName[3] = "KNN model - n. of positives";
        this.ResultsName[4] = "KNN model - n. of negatives";
        this.ResultsName[5] = "GBT model predicted activity";
        this.ResultsName[6] = "GBT model - confidence";
        this.ResultsName[7] = "Sarpy model predicted activity";

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
            int[] cons_res = new int[3];

            GBTResult res_gbt = model_gbt.predict(CurMolecule.getInputSMILES());
            cons_res[0] = res_gbt.prediction == Descriptor.MISSING_VALUE ? -1 : res_gbt.prediction;

            DioMolecularKNNClassifierFixed.Prediction res_knn = model_knn.predict(CurMolecule.getInputSMILES());
            try {
                cons_res[1] = Integer.parseInt(res_knn.prediction); // should be "1", "0", or "out of domain"
            } catch (Exception e) {
                cons_res[1] = -1;
            }

            int res_sarpy = model_sarpy.match(CurMolecule);
            cons_res[2] = res_sarpy;

            // majority votes
            int pos = 0, neg = 0;
            for (int pred : cons_res) {
                if (pred == 1)
                    pos++;
                else if (pred == 0)
                    neg++;
            }

            int final_res;

            // Majority votes
//            if (pos > neg)
//                final_res = 1;
//            else if (neg > pos)
//                final_res = 0;
//            else
//                final_res = -1;

            // Alternate scheme
            if (pos > 0)
                final_res = 1;
            else if (neg > 0)
                final_res = 0;
            else
                final_res = -1;


            CurOutput.setMainResultValue(final_res);

            String[] Res = new String[ResultsSize];
            try {
                Res[0] = this.GetTrainingSet().getClassLabel(final_res);
            } catch (Throwable ex) {
                log.warn("Unable to find label for DIO1 value " + final_res);
                Res[0] = Double.toString(final_res);
            }
            try {
                Res[1] = this.GetTrainingSet().getClassLabel(cons_res[1]);
            } catch (Throwable ex) {
                log.warn("Unable to find label for DIO1 value " + cons_res[1]);
                Res[1] = Double.toString(cons_res[1]);
            }
            Res[2] = "" + res_knn.similarity;
            Res[3] = "" + res_knn.positives;
            Res[4] = "" + res_knn.negatives;
            try {
                Res[5] = this.GetTrainingSet().getClassLabel(cons_res[0]);
            } catch (Throwable ex) {
                log.warn("Unable to find label for DIO1 value " + cons_res[0]);
                Res[5] = Double.toString(cons_res[0]);
            }
            Res[6] = "" + res_gbt.confidence;
            try {
                Res[7] = this.GetTrainingSet().getClassLabel(res_sarpy);
            } catch (Throwable ex) {
                log.warn("Unable to find label for DIO1 value " + res_sarpy);
                Res[7] = Double.toString(res_sarpy);
            }

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
