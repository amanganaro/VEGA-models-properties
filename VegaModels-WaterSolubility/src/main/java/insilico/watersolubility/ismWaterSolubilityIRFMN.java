package insilico.watersolubility;

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
import insilico.watersolubility.descriptors.EmbeddedDescriptors;
import lombok.extern.log4j.Log4j;
//import lombok.extern.log4j.Log4j;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

@Log4j
public class ismWaterSolubilityIRFMN extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_watersolubility_irfmn.xml";

    private static final double[] Nnet_norm_center = {
            2.10838614848032,
            3.26609591429995,
            1.3461649227703,
            -0.001482311908321,
            0.027355007473842,
            19.3401679123069,
            2.0865615346288,
            0.194070752366717,
            1.01420029895366,
            2.95398928749377,
            0.801831091180867,
            0.790287244643747,
            1.21076233183856,
            0.580312655705032,
            0.410812157448929
    };

    private static final double[] Nnet_norm_scale = {
            2.1030866739873,
            2.65146339266069,
            1.06371686135532,
            0.408898752071144,
            0.016914415799444,
            27.4034698946927,
            2.55948327870098,
            0.642737976651477,
            1.73184870037676,
            0.558366924041694,
            0.166154443659281,
            0.325812692726696,
            1.51409155232518,
            5.21746373654549,
            1.218837398151
    };

    private ModelANNFromPMML ANN;
    private double MW;

    public ismWaterSolubilityIRFMN()
            throws InitFailureException {
        super(ModelData);

        // Init PMML ANN
        try {
            URL src = getClass().getResource("/data/ws_nnet.pmml");
            ANN = new ModelANNFromPMML(src.openStream(), "WS.value");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }


        // Define no. of descriptors
        this.DescriptorsSize = 15;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ALogP";
        this.DescriptorsNames[1] = "piPC8";
        this.DescriptorsNames[2] = "SEigm";
        this.DescriptorsNames[3] = "MATS1v";
        this.DescriptorsNames[4] = "PW6";
        this.DescriptorsNames[5] = "P_VSA_i_4";
        this.DescriptorsNames[6] = "LogP";
        this.DescriptorsNames[7] = "nR10";
        this.DescriptorsNames[8] = "CATS2D_4_DL";
        this.DescriptorsNames[9] = "BEH2p";
        this.DescriptorsNames[10] = "BIC3";
        this.DescriptorsNames[11] = "GATS2m";
        this.DescriptorsNames[12] = "H.050";
        this.DescriptorsNames[13] = "D.Dr3";
        this.DescriptorsNames[14] = "CATS2D_7_DL";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted water solubility [-log(mol/L)]";
        this.ResultsName[1] = "Predicted water solubility [mg/L]";
        this.ResultsName[2] = "Molecular Weight";
        this.ResultsName[3] = "Experimental value [mg/L]";

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
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {

        try {
            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();
            embeddedDescriptors.CalculateAllEmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.ALogP;
            Descriptors[1] = embeddedDescriptors.piPC8;
            Descriptors[2] = embeddedDescriptors.SEigm;
            Descriptors[3] = embeddedDescriptors.MATS1v;
            Descriptors[4] = embeddedDescriptors.PW6;
            Descriptors[5] = embeddedDescriptors.P_VSA_I_4;
            Descriptors[6] = embeddedDescriptors.LogP;
            Descriptors[7] = embeddedDescriptors.nR10;
            Descriptors[8] = embeddedDescriptors.CATS2D_4_DL;
            Descriptors[9] = embeddedDescriptors.BEH2p;
            Descriptors[10] = embeddedDescriptors.BIC3;
            Descriptors[11] = embeddedDescriptors.GATS2m;
            Descriptors[12] = embeddedDescriptors.H050;
            Descriptors[13] = embeddedDescriptors.DDr3;
            Descriptors[14] = embeddedDescriptors.CATS2D_7_DL;

            MW = CurMolecule.GetBasicDescriptorByName("MW_da").getValue();
        } catch (Throwable e) {
            log.warn(e.getMessage());
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
            double RoundedVal = (double)Math.round(Descriptors[i] * 1000d) / 1000d;
            NormDescriptors[i] = (RoundedVal - Nnet_norm_center[i]) / Nnet_norm_scale[i];
        }

//            NormDescriptors[i] = (Descriptors[i] - Nnet_norm_center[i]) / Nnet_norm_scale[i];

        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        argumentsObject.put("ALogP", NormDescriptors[0]);
        argumentsObject.put("piPC8", NormDescriptors[1]);
        argumentsObject.put("SEigm", NormDescriptors[2]);
        argumentsObject.put("MATS1v", NormDescriptors[3]);
        argumentsObject.put("PW6", NormDescriptors[4]);
        argumentsObject.put("P_VSA_i_4", NormDescriptors[5]);
        argumentsObject.put("LogP", NormDescriptors[6]);
        argumentsObject.put("nR10", NormDescriptors[7]);
        argumentsObject.put("CATS2D_4_DL", NormDescriptors[8]);
        argumentsObject.put("BEH2p", NormDescriptors[9]);
        argumentsObject.put("BIC3", NormDescriptors[10]);
        argumentsObject.put("GATS2m", NormDescriptors[11]);
        argumentsObject.put("H.050", NormDescriptors[12]);
        argumentsObject.put("D.Dr3", NormDescriptors[13]);
        argumentsObject.put("CATS2D_7_DL", NormDescriptors[14]);

        double Prediction;
        try {
            Prediction = ANN.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // −Log10(mol/L)

        double ConvertedValue = Math.pow(10, -Prediction) * MW * 1000;
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // in mg/L
        else if (ConvertedValue>0.0001)
            Res[1] = Format_4D.format(ConvertedValue); // in mg/L
        else
            Res[1] = Format_6D.format(ConvertedValue); // in mg/L

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
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.2, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.2, 0.6);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.2, 0.6);
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
            double ConvertedValue = Math.pow(10, -CurOutput.getExperimental()) * MW * 1000;
            if (ConvertedValue>1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedValue); // in mg/L
            else if (ConvertedValue>0.0001)
                CurOutput.getResults()[3] = Format_4D.format(ConvertedValue); // in mg/L
            else
                CurOutput.getResults()[3] = Format_6D.format(ConvertedValue); // in mg/L
        }

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

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

        // Always gray light - no threshold for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }

}
