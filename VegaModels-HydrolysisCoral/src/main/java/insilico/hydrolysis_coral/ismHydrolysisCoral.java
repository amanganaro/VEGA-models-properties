package insilico.hydrolysis_coral;

import insilico.coral.CoralModel;
import insilico.coral.models.hydrolysis.CoralHydrolysis;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismHydrolysisCoral extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_hydrolysis_coral.xml";

    
    private CoralModel Hydro;
    
    
    public ismHydrolysisCoral() 
            throws InitFailureException {
        super(ModelData);
        
        // Init Coral model
        try {  
            Hydro = new CoralHydrolysis();
        } catch (Exception ex) {
            throw new InitFailureException("Unable to init coral model");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Hydrolysis [log(days)]";
        this.ResultsName[1] = "Predicted Hydrolysis [days]";
        this.ResultsName[2] = "Experimental Hydrolysis [days]";
        
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
        return blocks;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
                       
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {
        
        double Prediction;
        try {
            Prediction = Hydro.Predict(this.CurMolecule.GetSMILES());
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_3D.format(Prediction)); // log(days)
        double ConvertedValue = Math.pow(10, Prediction);
        Res[1] = Format_2D.format(ConvertedValue); // days
        Res[2] = "-";
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
        
        // Add transformed (mg/kg) experimental if needed
        if (CurOutput.HasExperimental()) {
            double ConvertedValue = Math.pow(10, CurOutput.getExperimental());
            CurOutput.getResults()[2] = Format_2D.format(ConvertedValue); // days
        }
             
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
       
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // No assessment color for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
