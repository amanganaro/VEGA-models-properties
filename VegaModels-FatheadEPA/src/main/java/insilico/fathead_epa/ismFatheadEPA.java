package insilico.fathead_epa;

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
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.fathead_epa.descriptors.EmbeddedDescriptors;
import insilico.fathead_epa.descriptors.FatheadEPARegression;
import insilico.fathead_epa.descriptors.LC50FMDescriptors;

import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismFatheadEPA extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_fathead_epa.xml";
    
    private final LC50FMDescriptors CustomDescriptors;
    private double MW;
    
    
    public ismFatheadEPA() 
            throws InitFailureException {
        super(ModelData);
        
        // Custom descriptors
        CustomDescriptors = new LC50FMDescriptors();
        
        // Define no. of descriptors
        this.DescriptorsSize = 21;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "SsssN";
        this.DescriptorsNames[1] = "SdssS_acnt";
        this.DescriptorsNames[2] = "MDE_C_33";
        this.DescriptorsNames[3] = "MDE_O_11";
        this.DescriptorsNames[4] = "BEH2m";
        this.DescriptorsNames[5] = "nDblBo";
        this.DescriptorsNames[6] = "nS";
        this.DescriptorsNames[7] = "nR9";
        this.DescriptorsNames[8] = "ATS5p";
        this.DescriptorsNames[9] = "MATS7e";
        this.DescriptorsNames[10] = "GATS3e";
        this.DescriptorsNames[11] = "SRW3";
        this.DescriptorsNames[12] = "ALogP";
        this.DescriptorsNames[13] = "n_aNO2";
        this.DescriptorsNames[14] = "n_CHO";
        this.DescriptorsNames[15] = "n_POH";
        this.DescriptorsNames[16] = "n_NNO";
        this.DescriptorsNames[17] = "n_NCON";
        this.DescriptorsNames[18] = "n_CN";
        this.DescriptorsNames[19] = "n_Cl";
        this.DescriptorsNames[20] = "n_aCHO";
        
        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted toxicity [-log(mol/l)]";
        this.ResultsName[1] = "Predicted toxicity [mg/l]";
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
            
            // Some custom descriptors
            CustomDescriptors.Calculate(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);


            Descriptors[0] = embeddedDescriptors.getSsssN();
            Descriptors[1] = CustomDescriptors.getN_sssS();
            Descriptors[2] = embeddedDescriptors.getMDE_C_33();
            Descriptors[3] = embeddedDescriptors.getMDE_O_11();
            Descriptors[4] = embeddedDescriptors.getBEH2m();
            Descriptors[5] = embeddedDescriptors.getNDblBo();
            Descriptors[6] = embeddedDescriptors.getNS();
            Descriptors[7] = embeddedDescriptors.getNR9();
            Descriptors[8] = embeddedDescriptors.getATS5p();
            Descriptors[9] = embeddedDescriptors.getMATS7e();
            Descriptors[10] = embeddedDescriptors.getGATS3e();
            Descriptors[11] = embeddedDescriptors.getSRW3();
            Descriptors[12] = embeddedDescriptors.getALogP();
            Descriptors[13] = CustomDescriptors.getN_aNO2();
            Descriptors[14] = CustomDescriptors.getN_CHO();
            Descriptors[15] = CustomDescriptors.getN_POH();
            Descriptors[16] = CustomDescriptors.getN_NNO();
            Descriptors[17] = CustomDescriptors.getN_NCON();
            Descriptors[18] = CustomDescriptors.getN_CN();
            Descriptors[19] = CustomDescriptors.getN_Cl();
            Descriptors[20] = CustomDescriptors.getN_aCHO();

            MW = CurMolecule.GetBasicDescriptorByName("MW_da").getValue();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        double Prediction = FatheadEPARegression.Calculate(Descriptors);
        
        CurOutput.setMainResultValue(Prediction);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // toxicity -log(mol/L)
        double ConvertedValue = Math.pow(10, (-1 * Prediction)) * 1000 * MW;
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // in mg/L
        else
            Res[1] = Format_4D.format(ConvertedValue); // in mg/L
        Res[2] = Format_2D.format(MW); // MW
        Res[3] = "-"; // Converted experimental - set after
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
            double ConvertedExp = Math.pow(10, (-1 * CurOutput.getExperimental())) * 1000 * MW;
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
        Val = Math.pow(10, (-1 * Val)) * 1000 * MW; // convert to mg/L
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
