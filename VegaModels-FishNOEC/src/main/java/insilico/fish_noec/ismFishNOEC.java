package insilico.fish_noec;

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
import insilico.fish_noec.descriptors.EmbeddedDescriptors;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
@Slf4j
public class ismFishNOEC extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_fish_noec.xml";
    
    private static final double[] NormalizationStd = {1.23886463598987,1.6112255765033,2.95431218043379,2.51194581804481,0.973281400564443};
    private static final double[] NormalizationMean = {1.71518666666667,3.09192,4.8548,4.84473333333333,2.28166666666667};

    
    private ModelANNFromPMML Model;
    private double MW;
    
    
    public ismFishNOEC() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML model
        try {  
            URL src = getClass().getResource("/data/fish_noec_model.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 5;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ATS5p";
        this.DescriptorsNames[1] = "ALogP";
        this.DescriptorsNames[2] = "X3sol";
        this.DescriptorsNames[3] = "X1v";
        this.DescriptorsNames[4] = "BEH4m";
        
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
            
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors[0] = embeddedDescriptors.getATS5p();
            Descriptors[1] = embeddedDescriptors.getALogP();
            Descriptors[2] = embeddedDescriptors.getX3sol();
            Descriptors[3] = embeddedDescriptors.getX1v();
            Descriptors[4] = embeddedDescriptors.getBEH4m();
            
            // MW in constitutional is given as a SCALED 
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            MW = embeddedDescriptors.getMW() * CarbonWeight;
            
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
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // box-cox, lambda = 0.01
        double ConvertedValue = Math.pow( (Prediction * 0.01 + 1), (1.0 / 0.01) ) * MW;
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
            double ConvertedValue = Math.pow( (CurOutput.getExperimental() * 0.01 + 1), (1.0 / 0.01) ) * MW;
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
        Val = Math.pow( (Val * 0.01 + 1), (1.0 / 0.01) ) * MW;
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
