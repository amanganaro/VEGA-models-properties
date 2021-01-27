package insilico.carcinogenicity_sfiregression;

import insilico.carcinogenicity_sfiregression.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
@Slf4j
public class ismCarcinogenicitySFIRegression extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_carc_sfiregression.xml";
    
    private ModelANNFromPMML ANN;
    
    
    public ismCarcinogenicitySFIRegression() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML ANN
        try {  
            URL src = getClass().getResource("/data/sfi_ann.pmml");
            ANN = new ModelANNFromPMML(src.openStream(), "Log SF");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 12;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "C-041";
        this.DescriptorsNames[1] = "ATS8m";
        this.DescriptorsNames[2] = "GATS6p";
        this.DescriptorsNames[3] = "CATS2D_3_DL";
        this.DescriptorsNames[4] = "CATS2D_7_DL";
        this.DescriptorsNames[5] = "nN-N";
        this.DescriptorsNames[6] = "IC4";
        this.DescriptorsNames[7] = "B2(Cl..Cl)";
        this.DescriptorsNames[8] = "B4(O..Cl)";
        this.DescriptorsNames[9] = "B7(Cl..Cl)";
        this.DescriptorsNames[10] = "B8(Cl..Cl)";
        this.DescriptorsNames[11] = "SRW7";
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Inhalation Carcinogenicity SF (log form) [log(1/(mg/kg-day))]";
        this.ResultsName[1] = "Predicted Inhalation Carcinogenicity SF [1/(mg/kg-day)]";
        this.ResultsName[2] = "Experimental value [1/(mg/kg-day)]";
        
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

            Descriptors[0] = embeddedDescriptors.C_041;
            Descriptors[1] = embeddedDescriptors.ATS8m;
            Descriptors[2] = embeddedDescriptors.GATS6p;
            Descriptors[3] = embeddedDescriptors.CATS2D_3_DL;
            Descriptors[4] = embeddedDescriptors.CATS2D_7_DL;
            Descriptors[5] = embeddedDescriptors.nN_N;
            Descriptors[6] = embeddedDescriptors.IC4;
            Descriptors[7] = embeddedDescriptors.B2Cl_Cl;
            Descriptors[8] = embeddedDescriptors.B4O_Cl;
            Descriptors[9] = embeddedDescriptors.B7Cl_Cl;
            Descriptors[10] = embeddedDescriptors.B8Cl_Cl;
            Descriptors[11] = embeddedDescriptors.SRW7;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        argumentsObject.put("C-041", Descriptors[0]);
        argumentsObject.put("ATS8m", Descriptors[1]);
        argumentsObject.put("GATS6p", Descriptors[2]);
        argumentsObject.put("CATS2D_3_DL", Descriptors[3]);
        argumentsObject.put("CATS2D_7_DL", Descriptors[4]);
        argumentsObject.put("nN-N", Descriptors[5]);
        argumentsObject.put("IC4", Descriptors[6]);
        argumentsObject.put("B2(Cl..Cl)", Descriptors[7]);
        argumentsObject.put("B4(O..Cl)", Descriptors[8]);
        argumentsObject.put("B7(Cl..Cl)", Descriptors[9]);
        argumentsObject.put("B8(Cl..Cl)", Descriptors[10]);
        argumentsObject.put("SRW7", Descriptors[11]);
        
        double PredictionNotNorm;
        try {
            PredictionNotNorm = ANN.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        // Result is normalized with the params used in the original ANN
        double scale = 0.1104129960968328;
        double translate = 0.49199279762690556;
        double Prediction = (1.0 / scale) * PredictionNotNorm - (translate / scale);
        
        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // log(1/(mg/kg-day))
        double ConvertedValue = Math.pow(10, Prediction);
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // 1/(mg/kg-day)
        else
            Res[1] = Format_4D.format(ConvertedValue); // 1/(mg/kg-day)
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
            double ConvertedValue = Math.pow(10, CurOutput.getExperimental());
            if (ConvertedValue>1)
                CurOutput.getResults()[2] = Format_2D.format(ConvertedValue); // 1/(mg/kg-day)
            else
                CurOutput.getResults()[2] = Format_4D.format(ConvertedValue); // 1/(mg/kg-day)
        }
                
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
       
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Always gray light - no threshold for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
