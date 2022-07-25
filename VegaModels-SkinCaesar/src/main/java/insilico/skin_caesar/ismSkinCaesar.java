package insilico.skin_caesar;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.skin_caesar.descriptors.EmbeddedDescriptors;
import lombok.extern.log4j.Log4j;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class ismSkinCaesar extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_skin_caesar.xml";
    
    private final ModelAFP AFP;
    private ModelAFPResult AFPResults;
    
    
    public ismSkinCaesar() 
            throws InitFailureException {
        super(ModelData);
        
        // Builds model objects
        AFP = new ModelAFP("/data/skin_AFP_rules.csv");
        
        // Define no. of descriptors
        this.DescriptorsSize = 7;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nN";
        this.DescriptorsNames[1] = "GNar";
        this.DescriptorsNames[2] = "X2v";
        this.DescriptorsNames[3] = "EEig10ri";
        this.DescriptorsNames[4] = "GGI8";
        this.DescriptorsNames[5] = "nCconj";
        this.DescriptorsNames[6] = "O-058";
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted skin sensitization activity";
        this.ResultsName[1] = "O(Active)";
        this.ResultsName[2] = "O(Inactive)";
        
        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);


            Descriptors[0] = embeddedDescriptors.getNN();
            Descriptors[1] = embeddedDescriptors.getGNar();
            Descriptors[2] = embeddedDescriptors.getX2v();
            Descriptors[3] = embeddedDescriptors.getEEig10ri();
            Descriptors[4] = embeddedDescriptors.getGGI8();
            Descriptors[5] = embeddedDescriptors.getNCconj();
            Descriptors[6] = embeddedDescriptors.getO_058();
                        
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        AFPResults = AFP.Calculate_Prediction(Descriptors);
        
        int MainResult;
        
        if (AFPResults.ResultVal == 2) {
            MainResult = 1;
        } else if (AFPResults.ResultVal == 1) {
            MainResult = 0;
        } else {
            MainResult = -1;
        }
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // skin classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for skin value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Format_2D.format(AFPResults.O_Active); // O value for active
        Res[2] = Format_2D.format(AFPResults.O_Inactive); // O value for inactive
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(2);

        // (only retrieve similar molecules if n.a. prediction)
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == -1) {
            try {
                adq.SetSimilarMolecules(CurMolecule, CurOutput);
            } catch (GenericFailureException ex) {
                // do nothing
            }
            return InsilicoModel.AD_ERROR;
        }

        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.8, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.8, 0.6);
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

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.8, 0.6);
        CurOutput.setADI(ADI);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        
    }
    
}
