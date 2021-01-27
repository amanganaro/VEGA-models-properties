package insilico.algae_ec50;

import insilico.algae_ec50.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
//import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

public class ismAlgaeEC50 extends InsilicoModel {
    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_algae_ec50.xml";

    private static final double[] NormalizationStd = {1.845303661, 1.399946997, 1.013529117, 0.409783504, 0.99737062, 0.332882058, 1.009054835, 5.129926891, 1.535201118, 34.81612088};
    private static final double[] NormalizationMean = {1.920075397, 2.074214286, 1.068178571, -0.199599206, 0.839785714, -0.136900794, 1.713869048, 13.08706349, 0.924603175, 16.23630556};


    private ModelANNFromPMML Model;
    private double MW;


    public ismAlgaeEC50()
            throws InitFailureException {
        super(ModelData);

        // Init PMML model
        try {
            URL src = getClass().getResource("/data/algae_ec50_model.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 10;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ALogP";
        this.DescriptorsNames[1] = "MLogP";
        this.DescriptorsNames[2] = "BEH5p";
        this.DescriptorsNames[3] = "MATS4s";
        this.DescriptorsNames[4] = "BEL5e";
        this.DescriptorsNames[5] = "MATS3s";
        this.DescriptorsNames[6] = "X3v";
        this.DescriptorsNames[7] = "MW";
        this.DescriptorsNames[8] = "CATS2D_3_DL";
        this.DescriptorsNames[9] = "P_VSA_m_4";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted EC50 [a-dimensional]";
        this.ResultsName[1] = "Predicted EC50 [mg/l]";
        this.ResultsName[2] = "Molecular Weight";
        this.ResultsName[3] = "Experimental value [mg/l]";

        // Define AD items
        this.ADItemsName = new String[6];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexRange().GetIndexName();
        this.ADItemsName[5] = new ADIndexACF().GetIndexName();

    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();

            embeddedDescriptors.CalculateDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.ALogP;
            Descriptors[1] = embeddedDescriptors.MLogP;
            Descriptors[2] = embeddedDescriptors.BEH5p;
            Descriptors[3] = embeddedDescriptors.MATS4s;
            Descriptors[4] = embeddedDescriptors.BEL5e;
            Descriptors[5] = embeddedDescriptors.MATS3s;
            Descriptors[6] = embeddedDescriptors.X3v;
            Descriptors[7] = embeddedDescriptors.Mw;
            Descriptors[8] = embeddedDescriptors.CATS2D_3_DL;
            Descriptors[9] = embeddedDescriptors.P_VSA_m_4;

            // MW in constitutional is given as a SCALED
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            MW = CarbonWeight * embeddedDescriptors.Mw;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }


    @Override
    protected short CalculateModel() {

        // Scale input descriptors
        double[] ScaledDescriptors = new double[this.DescriptorsSize];
        for (int i=0; i<DescriptorsSize; i++)
            ScaledDescriptors[i] = (Descriptors[i] - NormalizationMean[i]) / NormalizationStd[i];

        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        for (int i=0; i<DescriptorsSize; i++)
            argumentsObject.put(this.DescriptorsNames[i], ScaledDescriptors[i]);

        // Run pmml model
        double Prediction;
        try {
            Prediction = Model.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // box-cox, lambda = 0.07
        double ConvertedValue = Math.pow( (Prediction * 0.07 + 1), (1.0 / 0.07) ) * MW;
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // mg/L
        else
            Res[1] = Format_4D.format(ConvertedValue); // mg/L
        Res[2] = Format_2D.format(MW); // MW
        Res[3] = "-";
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.5, 0.8);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.5, 0.8);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.5, 0.8);
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

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);

        // Add transformed (mg/L) experimental if needed
        if (CurOutput.HasExperimental()) {
            double ConvertedValue = Math.pow( (CurOutput.getExperimental() * 0.07 + 1), (1.0 / 0.07) ) * MW;
            if (ConvertedValue>1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedValue); // mg/l
            else
                CurOutput.getResults()[3] = Format_4D.format(ConvertedValue); // mg/l
        }

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {
        double LC_threshold_red = 1; // in mg/l
        double LC_threshold_orange = 10; // in mg/l
        double LC_threshold_yellow = 100; // in mg/l

        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set (mg/L) if available

        String ADItemWarnings =
                ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());

        String Result = CurOutput.getResults()[1] + " mg/L";

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
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L"));
        }

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        Val = Math.pow( (Val * 0.07 + 1), (1.0 / 0.07) ) * MW; // convert to mg/L
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
