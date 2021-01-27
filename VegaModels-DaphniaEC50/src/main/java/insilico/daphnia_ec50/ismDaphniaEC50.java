package insilico.daphnia_ec50;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import insilico.daphnia_ec50.descriptors.EmbeddedDescriptors;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
public class ismDaphniaEC50 extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_daphnia_ec50.xml";

    private double MW;

    private static final double[] Nnet_norm_center = {
        -0.48788141025641,
        2.33004807692308,
        0.538461538461538,
        0.025641025641026,
        -0.088173076923077,
        0.092948717948718,
        0.695810897435897,
        -0.106516025641026,
        0.019230769230769,
        0.061926282051282,
        1.00839423076923,
        0.019230769230769,
        -1.53240064102564
    };

    private static final double[] Nnet_norm_scale = {
        0.950344947718096,
        1.61924531097126,
        0.499319348734632,
        0.158315812370505,
        0.35147406930669,
        0.600922354234128,
        0.382159727565354,
        0.512425472331462,
        0.137555780963858,
        0.729688351977207,
        0.04385251918964,
        0.137555780963858,
        1.09135703077268
    };


    private ModelANNFromPMML ANN;


    public ismDaphniaEC50()
            throws InitFailureException {
        super(ModelData);

        // Init PMML ANN
        try {
            URL src = getClass().getResource("/data/daphnia_ec50_ann.pmml");
            ANN = new ModelANNFromPMML(src.openStream(), "log.mmo.l.");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 12;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "EEig15bo";
        this.DescriptorsNames[1] = "ALogP";
        this.DescriptorsNames[2] = "B2.C..O.";
        this.DescriptorsNames[3] = "S.106";
        this.DescriptorsNames[4] = "EEig8dm";
        this.DescriptorsNames[5] = "F4.Cl..Cl.";
        this.DescriptorsNames[6] = "GATS1m";
        this.DescriptorsNames[7] = "MATS4p";
        this.DescriptorsNames[8] = "B10.C..N.";
        this.DescriptorsNames[9] = "MATS5e";
        this.DescriptorsNames[10] = "Me";
        this.DescriptorsNames[11] = "F10.O..O.";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted daphnia EC50 (log form) [log(mmol/l))]";
        this.ResultsName[1] = "Predicted daphnia EC50 [mg/l]";
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

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.EEig15bo;
            Descriptors[1] = embeddedDescriptors.ALogP;
            Descriptors[2] = embeddedDescriptors.B2_C_O;
            Descriptors[3] = embeddedDescriptors.S_106;
            Descriptors[4] = embeddedDescriptors.EEig8dm;
            Descriptors[5] = embeddedDescriptors.F4_Cl_Cl;
            Descriptors[6] = embeddedDescriptors.GATS1m;
            Descriptors[7] = embeddedDescriptors.MATS4p;
            Descriptors[8] = embeddedDescriptors.B10_C_N;
            Descriptors[9] = embeddedDescriptors.MATS5e;
            Descriptors[10] = embeddedDescriptors.Me;
            Descriptors[11] = embeddedDescriptors.F10_O_O;

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

        // normalize descriptors for nnet
        double[] NormDescriptors = new double[DescriptorsSize];

        for (int i=0; i<DescriptorsSize; i++) {
            // rounded to 3rd dec digit to be consistent with descriptors and scaling in R
            double RoundedVal = (double) Math.round(Descriptors[i] * 1000d) / 1000d;
            NormDescriptors[i] = (RoundedVal - Nnet_norm_center[i]) / Nnet_norm_scale[i];
        }

        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        for (int i=0; i<DescriptorsSize; i++)
            argumentsObject.put(DescriptorsNames[i], NormDescriptors[i]);

        double Prediction;
        try {
            Prediction = ANN.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // log(mmol/l)
        double ConvertedValue = Math.pow(10, Prediction) * MW;
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // in mg/L
        else
            Res[1] = Format_4D.format(ConvertedValue); // in mg/L
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
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.2, 0.8);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.2, 0.8);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.2, 0.8);
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
            double ConvertedExp = Math.pow(10, CurOutput.getExperimental()) * MW;
            if (ConvertedExp>1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedExp); // in mg/L
            else
                CurOutput.getResults()[3] = Format_4D.format(ConvertedExp); // in mg/L
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
        Val = Math.pow(10, Val) * MW; // convert to mg/L
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
