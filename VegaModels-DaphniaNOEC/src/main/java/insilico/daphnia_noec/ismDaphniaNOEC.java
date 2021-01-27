package insilico.daphnia_noec;

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
import insilico.daphnia_noec.descriptors.EmbeddedDescriptors;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
public class ismDaphniaNOEC extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_daphnia_noec.xml";

    private double MW;
    
    private static final double[] Nnet_norm_center = {
        5.67037209302326,
        2.2555023255814,
        0.130232558139535,
        0.158139534883721,
        -0.039,
        0.02893023255814,
        -0.079548837209302,
        0.283720930232558,
        0.734883720930233,
        0.042944186046512
    };
    
    private static final double[] Nnet_norm_scale = {
        7.59548209732643,
        1.60154807956639,
        0.412160301203652,
        0.45663977814147,
        0.600114428340799,
        0.500362587938811,
        0.344815860632844,
        1.11425140046508,
        1.27500846999756,
        0.024230139784655        
    };    
    
    
    private ModelANNFromPMML ANN;
    
    
    public ismDaphniaNOEC() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML ANN
        try {  
            URL src = getClass().getResource("/data/daphnia_noec_ann.pmml");
            ANN = new ModelANNFromPMML(src.openStream(), "log.mmol.l.");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 10;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "OPerc";
        this.DescriptorsNames[1] = "MLogP";
        this.DescriptorsNames[2] = "nArNH2";
        this.DescriptorsNames[3] = "nS";
        this.DescriptorsNames[4] = "MATS5s";
        this.DescriptorsNames[5] = "MATS6m";
        this.DescriptorsNames[6] = "EEig7dm";
        this.DescriptorsNames[7] = "CATS2D_7_DL";
        this.DescriptorsNames[8] = "C.026";
        this.DescriptorsNames[9] = "JGI3";
        
        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted daphnia NOEC (log form) [log(mmol/l))]";
        this.ResultsName[1] = "Predicted daphnia NOEC [mg/l]";
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
            
            Descriptors[0] = embeddedDescriptors.OPerc;
            Descriptors[1] = embeddedDescriptors.MLogP;
            Descriptors[2] = embeddedDescriptors.nArNH2;
            Descriptors[3] = embeddedDescriptors.nS;
            Descriptors[4] = embeddedDescriptors.MATS5s;
            Descriptors[5] = embeddedDescriptors.MATS6m;
            Descriptors[6] = embeddedDescriptors.EEig7dm;
            Descriptors[7] = embeddedDescriptors.CATS2D_7_DL;
            Descriptors[8] = embeddedDescriptors.C_026;
            Descriptors[9] = embeddedDescriptors.JGI3;
            
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
        argumentsObject.put("OPerc", NormDescriptors[0]);
        argumentsObject.put("MLogP", NormDescriptors[1]);
        argumentsObject.put("nArNH2", NormDescriptors[2]);
        argumentsObject.put("nS", NormDescriptors[3]);
        argumentsObject.put("MATS5s", NormDescriptors[4]);
        argumentsObject.put("MATS6m", NormDescriptors[5]);
        argumentsObject.put("EEig7dm", NormDescriptors[6]);
        argumentsObject.put("CATS2D_7_DL", NormDescriptors[7]);
        argumentsObject.put("C.026", NormDescriptors[8]);
        argumentsObject.put("JGI3", NormDescriptors[9]);
        
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
