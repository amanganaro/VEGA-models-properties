package insilico.ppara_up;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import insilico.ppara_up.descriptors.EmbeddedDescriptors;
import lombok.SneakyThrows;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dmg.pmml.FieldName;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;


public class ismPPARAUp extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismPPARAUp.class);

    private static final double[] mean = {3.2228, 0.163, 0.106};
    private static final double[] stdDeviation = {1.6031, 0.031, 0.4287};


    private static final String ModelData = "/data/model_pparaup.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;

    public ismPPARAUp() throws InitFailureException {
        super(ModelData);

        try {
            URL src = getClass().getResource("/data/RF_model_pparaup.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Predicted Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 3;
        this.DescriptorsNames = new String[DescriptorsSize];

        this.DescriptorsNames[0] = "TI2_L";
        this.DescriptorsNames[1] = "PW4";
        this.DescriptorsNames[2] = "CATS2D_07_NL";

        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted activity";
        this.ResultsName[1] = "Inactivity class Probability";
        this.ResultsName[2] = "Activity class Probability";

    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule, false);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.TI2_L;
            Descriptors[1] = embeddedDescriptors.PW4;
            Descriptors[2] = embeddedDescriptors.CATS2D_07_NL;

        } catch (Throwable e){
            log.warn(e.getMessage());
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @SneakyThrows
    @Override
    protected short CalculateModel() {
        Map<String, Object> argumentsObject = new LinkedHashMap<>();
        double[] ScaledDescriptors = new double[this.DescriptorsSize];

        for (int i=0; i<DescriptorsSize; i++) {
            ScaledDescriptors[i] = (Descriptors[i] - mean[i]) / stdDeviation[i];
        }
        try {
            for (int i=0; i<DescriptorsSize; i++)
                argumentsObject.put(this.DescriptorsNames[i], ScaledDescriptors[i]);
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

        // Run pmml model
        int Prediction;
        double prediction_0;
        double prediction_1;

        try {
//            Prediction = Model.Evaluate(argumentsObject);
            Map<FieldName, ?> outputs = Model.EvaluateFullOutput(argumentsObject);
            Prediction = Integer.parseInt((String) outputs.get(FieldName.create("Predicted_Exp")));
            prediction_0 = (Double) outputs.get(FieldName.create("Probability_0"));
            prediction_1 = (Double) outputs.get(FieldName.create("Probability_1"));
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction);
        } catch (Throwable ex) {
            log.warn("Unable to find label for value " + Prediction);
            Res[0] = Integer.toString(Prediction);
        }
        Res[1] = Format_3D.format(prediction_0);
        Res[2] = Format_3D.format(prediction_1);

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(2);

        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.8, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.8, 0.6);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets Range check
        ADCheckDescriptorRange adrc = new ADCheckDescriptorRange();
        if (!adrc.Calculate(TS, Descriptors, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.8, 0.6);
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
