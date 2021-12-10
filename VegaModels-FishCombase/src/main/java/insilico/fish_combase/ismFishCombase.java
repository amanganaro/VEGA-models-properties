package insilico.fish_combase;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.coral.CoralModel;
import insilico.core.coral.models.biocides.fish.OptimizedCoralFish;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismFishCombase extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_fish_combase.xml";
    
    
    private CoralModel FishModel;
    private double MW;
    
    
    public ismFishCombase() 
            throws InitFailureException {
        super(ModelData);
        
        // Init Coral model
        try {  
            FishModel = new OptimizedCoralFish();
        } catch (Exception ex) {
            throw new InitFailureException("Unable to init coral model");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LC50 [-log(mmol/L)]";
        this.ResultsName[1] = "Predicted LC50 [mg/l]";
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
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
        
        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        DescriptorBlock desc;
        
        // MW
        desc = new Constitutional();
        blocks.add(desc);

        return blocks;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
                       
            // MW in constitutional is given as a SCALED 
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            MW = CarbonWeight * CurMolecule.GetBasicDescriptorByName("MW").getValue();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        double Prediction;
        try {
            Prediction = FishModel.Predict(this.CurMolecule.GetSMILES());
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);


        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // -log(mmol/L)
        double ConvertedValue = Math.pow(10, -Prediction) * MW;
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
            double ConvertedValue = Math.pow(10, -CurOutput.getExperimental()) * MW;
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
        Val = Math.pow(10, -Val) * MW;
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
