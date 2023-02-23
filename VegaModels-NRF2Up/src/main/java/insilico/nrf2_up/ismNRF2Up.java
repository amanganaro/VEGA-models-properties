package insilico.nrf2_up;

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
import insilico.nrf2_up.descriptors.EmbeddedDescriptors;
import insilico.nrf2_up.descriptors.MLogP;
import lombok.SneakyThrows;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dmg.pmml.FieldName;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;


public class ismNRF2Up extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismNRF2Up.class);

    private static final double[] mean = {2.265985932,72.8069519343,83.3807702227,1.6320257913,77.8901148886,7.5287221571,2.4724501758,3.368042204,4.427901524,5.1652989449,4.0572450176,42.6454279015,62.3047221571};
    private static final double[] stdDeviation = {1.5379366535,41.9561232784,41.1852020728,0.232307716,38.9656361906,9.6637448823,2.1303206186,0.7763722853,5.8889185286,7.1158829198,0.6753593947,47.4645068958,41.8755590915};

    private boolean loadDescriptorsFromFile = false;

    private static final String ModelData = "/data/model_nrf2_up.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;


    public void setLoadDescriptorsFromFile(boolean fromFile){
        this.loadDescriptorsFromFile = fromFile;
    }

    public ismNRF2Up() throws InitFailureException {
        super(ModelData);

        try {
            URL src = getClass().getResource("/data/RF_model_nrf2up.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Predicted Exp");
        } catch (IOException | InitFailureException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 13;
        this.DescriptorsNames = new String[DescriptorsSize];

        this.DescriptorsNames[0] = "MLOGP";
        this.DescriptorsNames[1] = "P_VSA_ppp_L";
        this.DescriptorsNames[2] = "P_VSA_e_2";
        this.DescriptorsNames[3] = "SpMin2_Bh.e.";
        this.DescriptorsNames[4] = "P_VSA_i_2";
        this.DescriptorsNames[5] = "F07.C.C.";
        this.DescriptorsNames[6] = "NaasC";
        this.DescriptorsNames[7] = "IDDE";
        this.DescriptorsNames[8] = "CATS2D_06_LL";
        this.DescriptorsNames[9] = "CATS2D_04_LL";
        this.DescriptorsNames[10] = "IC3";
        this.DescriptorsNames[11] = "ATSC2s";
        this.DescriptorsNames[12] = "P_VSA_ppp_cyc";


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
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule, loadDescriptorsFromFile);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.MLOGP;
            Descriptors[1] = embeddedDescriptors.P_VSA_ppp_L;
            Descriptors[2] = embeddedDescriptors.P_VSA_e_2;
            Descriptors[3] = embeddedDescriptors.SpMin2_Bh_e;
            Descriptors[4] = embeddedDescriptors.P_VSA_i_2;
            Descriptors[5] = embeddedDescriptors.F07_C_C;
            Descriptors[6] = embeddedDescriptors.NaasC;
            Descriptors[7] = embeddedDescriptors.IDDE;
            Descriptors[8] = embeddedDescriptors.CATS2D_06_LL;
            Descriptors[9] = embeddedDescriptors.CATS2D_04_LL;
            Descriptors[10] = embeddedDescriptors.IC3;
            Descriptors[11] = embeddedDescriptors.ATSC2s;
            Descriptors[12] = embeddedDescriptors.P_VSA_ppp_cyc;

        } catch (Exception ex){
            log.warn(ex.getMessage());
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
