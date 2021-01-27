package insilico.daphnia_epa;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
//import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.daphnia_epa.descriptors.EmbeddedDescriptors;
import insilico.daphnia_epa.descriptors.LC50DMDescriptors;

import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismDaphniaEPA extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_daphnia_epa.xml";
    
    private final LC50DMDescriptors CustomDescriptors;
    private double MW;
    
    
    public ismDaphniaEPA() 
            throws InitFailureException {
        super(ModelData);
        
        // Custom descriptors
        CustomDescriptors = new LC50DMDescriptors();
        
        // Define no. of descriptors
        this.DescriptorsSize = 18;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "xc4";
        this.DescriptorsNames[1] = "StN";
        this.DescriptorsNames[2] = "SsSH";
        this.DescriptorsNames[3] = "SsOH_acnt";
        this.DescriptorsNames[4] = "HMax";
        this.DescriptorsNames[5] = "ssi";
        this.DescriptorsNames[6] = "MDE_N_33";
        this.DescriptorsNames[7] = "BEH1m";
        this.DescriptorsNames[8] = "BEH1p";
        this.DescriptorsNames[9] = "Mv";
        this.DescriptorsNames[10] = "MATS1m";
        this.DescriptorsNames[11] = "MATS1e";
        this.DescriptorsNames[12] = "GATS3m";
        this.DescriptorsNames[13] = "AMR";
        this.DescriptorsNames[14] = "CdS2N";
        this.DescriptorsNames[15] = "nAN";
        this.DescriptorsNames[16] = "nNP";
        this.DescriptorsNames[17] = "nSdOdO";
        
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
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = CustomDescriptors.getXc4();
            Descriptors[1] = embeddedDescriptors.StN;
            Descriptors[2] = embeddedDescriptors.SsSH;
            
            int SsOH_count = (int) embeddedDescriptors.nROH;
            SsOH_count += (int) embeddedDescriptors.nArOH;
            Descriptors[3] = SsOH_count;
            Descriptors[4] = embeddedDescriptors.HMax;

            // To be consistent with the original model, Hmax is set to 0.0
            // when it is not possibile to calculate it (compounds without 
            // hydrogen atoms)
            if (Descriptors[4] == Descriptor.MISSING_VALUE)
                Descriptors[4] = 0.0;
            
            // ssi not implemented yet! fix value (mean value in training set)
            Descriptors[5] = 0.85;
                    
            Descriptors[6] = embeddedDescriptors.MDE_N_33;
            Descriptors[7] = embeddedDescriptors.BEH1m;
            Descriptors[8] = embeddedDescriptors.BEH1p;
            Descriptors[9] = embeddedDescriptors.Mv;
            Descriptors[10] = embeddedDescriptors.MATS1m;
            Descriptors[11] = embeddedDescriptors.MATS1e;
            Descriptors[12] = embeddedDescriptors.GATS3m;
            Descriptors[13] = embeddedDescriptors.AMR;
            Descriptors[14] = CustomDescriptors.getN_CdS2N();
            Descriptors[15] = CustomDescriptors.getN_AN();
            Descriptors[16] = CustomDescriptors.getN_NP();
            Descriptors[17] = CustomDescriptors.getN_SdOdO();
            
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

        double Prediction = DaphniaEPARegression.Calculate(Descriptors);
        
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

        String ADItemWarnings = ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());
        
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
