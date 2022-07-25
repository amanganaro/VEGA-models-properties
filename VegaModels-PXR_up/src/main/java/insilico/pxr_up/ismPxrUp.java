package insilico.pxr_up;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import insilico.pxr_up.descriptors.EmbeddedDescriptors;
import lombok.SneakyThrows;
import lombok.extern.log4j.Log4j;
import org.dmg.pmml.FieldName;
import org.jpmml.evaluator.mining.MiningModelEvaluator;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

@Log4j
public class ismPxrUp extends InsilicoModel {

    private static final String ModelData = "/data/model_pxrup.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;

    private static final double[] normalizationCenter = {0.0299764454, 0.4776059957, 1.55192077087794, 2.4641905782, 8.06624411134904, 39.0678244111, 76.767329764454, -0.259103854389722, 4.0144571734, 7.10323768736617, 3.11351713062099, 1.81327194860814, 1.26654282655246, 12002.0278372591, 1.64783940042827};
    private static final double[] normalizationScale = {0.0302561321, 0.0294393304, 1.21377318129206, 1.3600946779, 8.17291955784308, 44.3464601839895, 46.5170499975922, 1.09192150904467, 0.769358465079047, 0.998769289424617, 0.472744171501964, 0.22172102080432, 0.0402121888990575, 63429.7241237166, 1.48126928129873};

    public ismPxrUp() throws InitFailureException {
        super(ModelData);

        try {
            URL src = getClass().getResource("/data/PXR_up_RF_model.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Predicted Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 15;
        this.DescriptorsNames = new String[DescriptorsSize];

        this.DescriptorsNames[0] = "ChiA_Dz.Z.";
        this.DescriptorsNames[1] = "Eta_alpha_A";
        this.DescriptorsNames[2] = "GGI3";
        this.DescriptorsNames[3] = "MAXDN";
        this.DescriptorsNames[4] = "MLOGP2";
        this.DescriptorsNames[5] = "P_VSA_ppp_D";
        this.DescriptorsNames[6] = "P_VSA_ppp_L";
        this.DescriptorsNames[7] = "SdssC";
        this.DescriptorsNames[8] = "SpMax_B.m.";
        this.DescriptorsNames[9] = "SpMax_B.s.";
        this.DescriptorsNames[10] = "SpMax4_Bh.m.";
        this.DescriptorsNames[11] = "SpMin2_Bh.m.";
        this.DescriptorsNames[12] = "SpPosA_B.v.";
        this.DescriptorsNames[13] = "Wap";
        this.DescriptorsNames[14] = "X5v";

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

            Descriptors[0] = embeddedDescriptors.ChiA_Dz_Z;
            Descriptors[1] = embeddedDescriptors.Eta_alpha_A;
            Descriptors[2] = embeddedDescriptors.GGI3;
            Descriptors[3] = embeddedDescriptors.MAXDN;
            Descriptors[4] = embeddedDescriptors.MLOGP2;
            Descriptors[5] = embeddedDescriptors.P_VSA_ppp_D;
            Descriptors[6] = embeddedDescriptors.P_VSA_ppp_L;
            Descriptors[7] = embeddedDescriptors.SdssC;
            Descriptors[8] = embeddedDescriptors.SpMax_B_m;
            Descriptors[9] = embeddedDescriptors.SpMax_B_s;
            Descriptors[10] = embeddedDescriptors.SpMax4_Bh_m;
            Descriptors[11] = embeddedDescriptors.SpMin2_Bh_m;
            Descriptors[12] = embeddedDescriptors.SpPosA_B_v;
            Descriptors[13] = embeddedDescriptors.Wap;
            Descriptors[14] = embeddedDescriptors.X5v;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @SneakyThrows
    @Override
    protected short CalculateModel() {

        double[] ScaledDescriptors = new double[this.DescriptorsSize];
        for (int i=0; i<DescriptorsSize; i++)
            ScaledDescriptors[i] = (Descriptors[i] - normalizationCenter[i]) / normalizationScale[i];

        Map<String, Object> argumentsObject = new LinkedHashMap<>();
        try {
            for (int i=0; i<DescriptorsSize; i++)
                argumentsObject.put(this.DescriptorsNames[i], ScaledDescriptors[i]);
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


        // Run pmml model
        int Prediction;
        double Probability0;
        double Probability1;
        try {
//            Prediction = Model.Evaluate(argumentsObject);
            Map<FieldName, ?> outputs = Model.EvaluateFullOutput(argumentsObject);
            Prediction = Integer.parseInt((String) outputs.get(FieldName.create("Predicted_Exp")));
            Probability0 = (Double) outputs.get(FieldName.create("Probability_0"));
            Probability1 = (Double) outputs.get(FieldName.create("Probability_1"));

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
        Res[1] = Format_3D.format(Probability0);
        Res[2] = Format_3D.format(Probability1);

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
