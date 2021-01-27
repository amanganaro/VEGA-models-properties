package insilico.persistence_quantative_water_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.cpannatnic.CPANNPredictor;
import insilico.cpannatnic.CPANNResults;
import insilico.persistence_quantative_water_irfmn.descriptors.EmbeddedDescriptors;

import java.net.URL;
import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismPersistenceWaterQuantitativeIrfmn extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_water_quantitative_irfmn.xml";
    
    private static final String CPANNWeights = "/data/water_modelweights.unw";
    private static final double[] DescScaling_Min = {0.156, 0, 0, -1, -1, 0, 0, 0, -2.503, 0, 0, 0, 0, 0, 0, -1.3};
    private static final double[] DescScaling_Max = {7.371, 27, 6, 1, 1, 2, 5.268, 2, 1.429, 4, 8, 1.29, 6.333, 46, 3, 3.401};
    

    
    private CPANNPredictor cpann;
    
    
    public ismPersistenceWaterQuantitativeIrfmn() 
            throws InitFailureException {
        super(ModelData);
        
        // Create and init CPANN object
        URL u = getClass().getResource(CPANNWeights);
        cpann = new CPANNPredictor();
        try {
            cpann.LoadCPANNObject(u.openStream());
        } catch (Exception ex) {
            throw new InitFailureException("Unable to load CPANN weights");
        }							
        
        // Define no. of descriptors
        this.DescriptorsSize = 16;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "MLogP";
        this.DescriptorsNames[1] = "F1(C..C)";
        this.DescriptorsNames[2] = "CATS2D_3_DL";
        this.DescriptorsNames[3] = "MATS8v";
        this.DescriptorsNames[4] = "MATS1e";
        this.DescriptorsNames[5] = "O.060";
        this.DescriptorsNames[6] = "MAXDP";
        this.DescriptorsNames[7] = "C.008";
        this.DescriptorsNames[8] = "MATS6p";
        this.DescriptorsNames[9] = "nR9";
        this.DescriptorsNames[10] = "F8(C..Cl)";
        this.DescriptorsNames[11] = "ATSC6e";
        this.DescriptorsNames[12] = "GATS6v";
        this.DescriptorsNames[13] = "F4(C..C)";
        this.DescriptorsNames[14] = "F3(N..Cl)";
        this.DescriptorsNames[15] = "MATS8m";

        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted HL [log(days)]";
        this.ResultsName[1] = "Predicted HL [days]";
        this.ResultsName[2] = "Experimental value [days]";
        
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
            
            Descriptors[0] = embeddedDescriptors.getMLogP();
            Descriptors[1] = embeddedDescriptors.getF1_C_C();
            Descriptors[2] = embeddedDescriptors.getCATS2D_3_DL();
            Descriptors[3] = embeddedDescriptors.getMATS8v();
            Descriptors[4] = embeddedDescriptors.getMATS1e();
            Descriptors[5] = embeddedDescriptors.getO_060();
            Descriptors[6] = embeddedDescriptors.getMAXDP();
            Descriptors[7] = embeddedDescriptors.getC_008();
            Descriptors[8] = embeddedDescriptors.getMATS6p();
            Descriptors[9] = embeddedDescriptors.getNR9();
            Descriptors[10] = embeddedDescriptors.getF8_C_Cl();
            Descriptors[11] = embeddedDescriptors.getATSC6e();
            Descriptors[12] = embeddedDescriptors.getGATS6v();
            Descriptors[13] = embeddedDescriptors.getF4_C_C();
            Descriptors[14] = embeddedDescriptors.getF3_N_Cl();
            Descriptors[15] = embeddedDescriptors.getMATS8m();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        // normalize descriptors
        double[] norm_descriptors = new double[DescriptorsSize];
        for (int i=0; i<DescriptorsSize; i++)
            norm_descriptors[i] = (Descriptors[i] - DescScaling_Min[i] ) / ( DescScaling_Max[i] - DescScaling_Min[i] );
        
        // predict with CPANN
        cpann.LoadSingleCompound(norm_descriptors, DescriptorsNames);
        ArrayList<CPANNResults> res = cpann.PerformPrediction();
        double Prediction = res.get(0).getPrediction();
        
        CurOutput.setMainResultValue(Prediction);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // hl (log days)
        double ConvertedValue = Math.pow(10, (Prediction));
        Res[1] = Long.toString( Math.round(ConvertedValue) ); // in days
        Res[2] = "-"; // Converted experimental - set after
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
            double ConvertedExp = Math.pow(10, (CurOutput.getExperimental()));
            CurOutput.getResults()[2] = Long.toString( Math.round(ConvertedExp) ); // in days
        }
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        double threshold_P = 40; 
        double threshold_vP = 60; 
        
        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set (days) if available

        String ADItemWarnings =
                ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());
        
        String Result = CurOutput.getResults()[1] + " days";
        
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
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getResults()[2] + " days", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getResults()[2] + " days"));
        }
        

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        Val = Math.pow(10, (Val)); // convert to days
        if (Val < threshold_P)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val < threshold_vP)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);

    }
}
