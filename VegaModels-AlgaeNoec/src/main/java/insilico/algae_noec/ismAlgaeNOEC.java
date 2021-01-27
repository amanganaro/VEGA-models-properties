package insilico.algae_noec;

import insilico.algae_noec.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.ADIndex;
import insilico.core.ad.item.ADIndexACF;
import insilico.core.ad.item.ADIndexADIAggregate;
import insilico.core.ad.item.ADIndexAccuracy;
import insilico.core.ad.item.ADIndexConcordance;
import insilico.core.ad.item.ADIndexMaxError;
import insilico.core.ad.item.ADIndexRange;
import insilico.core.ad.item.ADIndexSimilarity;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
public class ismAlgaeNOEC extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_algae_noec.xml";

    private static final double[] NormalizationStd = {1.96580628946963, 1.57270679007972, 0.407508563802331, 1.03131710701777, 2.07719955834384, 32.8951110912294, 2.53471963094721, 0.221410388880736, 0.420070597187593, 0.450166150332932, 0.074720585843248, 4.01795287765648, 0.421206580880607, 0.769309315072312, 0.209147840709661, 1.48747286370226, 0.17838451166724, 43.5600847116374, 17.7584167812485, 0.032927333533898, 1.80592869978331, 0.533378608054094};
    private static final double[] NormalizationMean = {2.05642682926829, 2.13852743902439, 0.031094512195122, 0.738807926829268, 4.08796951219512, 54.4145457317073, 1.54493597560976, -0.112451219512195, 0.740356707317073, 0.103844512195122, 1.22189329268293, 4.84451219512195, 1.00453048780488, 2.43940243902439, 1.41829268292683, 1.19512195121951, 0.691198170731707, 18.1890243902439, 11.7759329268293, 0.060329268292683, 3.08689329268293, -0.071612804878049};

    
    private ModelANNFromPMML Model;
    private double MW;
    
    
    public ismAlgaeNOEC() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML model
        try {  
            URL src = getClass().getResource("/data/algae_noec_model.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 22;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ALogP";
        this.DescriptorsNames[1] = "MLogP";
        this.DescriptorsNames[2] = "MATS1p";
        this.DescriptorsNames[3] = "BEH6p";
        this.DescriptorsNames[4] = "X1v";
        this.DescriptorsNames[5] = "P_VSA_i_2";
        this.DescriptorsNames[6] = "piPC10";
        this.DescriptorsNames[7] = "MATS1s";
        this.DescriptorsNames[8] = "GATS2v";
        this.DescriptorsNames[9] = "MATS2v";
        this.DescriptorsNames[10] = "SpPosA_v";
        this.DescriptorsNames[11] = "H.047";
        this.DescriptorsNames[12] = "GATS3s";
        this.DescriptorsNames[13] = "BEH3m";
        this.DescriptorsNames[14] = "SpMAD_m";
        this.DescriptorsNames[15] = "nCp";
        this.DescriptorsNames[16] = "BIC1";
        this.DescriptorsNames[17] = "T.C..N.";
        this.DescriptorsNames[18] = "ATSC3s";
        this.DescriptorsNames[19] = "PW5";
        this.DescriptorsNames[20] = "X2v";
        this.DescriptorsNames[21] = "MATS4v";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted NOEC [a-dimensional]";
        this.ResultsName[1] = "Predicted NOEC [mg/l]";
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

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];
            
            Descriptors[0] = embeddedDescriptors.ALogP;
            Descriptors[1] = embeddedDescriptors.MLogP;
            Descriptors[2] = embeddedDescriptors.MATS1p;
            Descriptors[3] = embeddedDescriptors.BEH6p;
            Descriptors[4] = embeddedDescriptors.X1v;
            Descriptors[5] = embeddedDescriptors.P_VSA_i_2;
            Descriptors[6] = embeddedDescriptors.piPC10;
            Descriptors[7] = embeddedDescriptors.MATS1s;
            Descriptors[8] = embeddedDescriptors.GATS2v;
            Descriptors[9] = embeddedDescriptors.MATS2v;
            Descriptors[10] = embeddedDescriptors.SpPosA_v;
            Descriptors[11] = embeddedDescriptors.H047;
            Descriptors[12] = embeddedDescriptors.GATS3s;
            Descriptors[13] = embeddedDescriptors.BEH3m;
            Descriptors[14] = embeddedDescriptors.SpMAD_m;
            Descriptors[15] = embeddedDescriptors.nCp;
            Descriptors[16] = embeddedDescriptors.BIC1;
            Descriptors[17] = embeddedDescriptors.TCN;
            Descriptors[18] = embeddedDescriptors.ATSC3s;
            Descriptors[19] = embeddedDescriptors.PW5;
            Descriptors[20] = embeddedDescriptors.X2v;
            Descriptors[21] = embeddedDescriptors.MATS4v;



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
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // box-cox, lambda = 0.03
        double ConvertedValue = Math.pow( (Prediction * 0.03 + 1), (1.0 / 0.03) ) * MW;
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
            double ConvertedValue = Math.pow( (CurOutput.getExperimental() * 0.03 + 1), (1.0 / 0.03) ) * MW;
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
        Val = Math.pow( (Val * 0.03 + 1), (1.0 / 0.03) ) * MW;  // convert to mg/L
        if (Val < LC_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val < LC_threshold_orange)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        else if (Val < LC_threshold_yellow)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);    }
}
