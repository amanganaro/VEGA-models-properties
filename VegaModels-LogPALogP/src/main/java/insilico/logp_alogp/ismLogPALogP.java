package insilico.logp_alogp;

import insilico.core.tools.utils.ModelUtilities;
import insilico.meylanlogp.descriptors.ADIndexADILogP;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.ADIndexADI;
import insilico.core.ad.item.ADIndexAccuracy;
import insilico.core.ad.item.ADIndexConcordance;
import insilico.core.ad.item.ADIndexMaxError;
import insilico.core.ad.item.ADIndexSimilarity;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.logp_alogp.descriptors.EmbeddedDescriptors;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismLogPALogP extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_logp_alogp.xml";
    
    private boolean SkipExperimentalInAD;
    
    public ismLogPALogP() 
            throws InitFailureException {
        super(ModelData);
        
        SkipExperimentalInAD = false;
        
        // Define no. of descriptors
        this.DescriptorsSize = 1;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ALogP";
        
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LogP";
        
        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        
    }
    
    public void SetSkipExperimentalInAD(boolean skip) {
        this.SkipExperimentalInAD = skip;
    }
    

    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];
            Descriptors[0] = embeddedDescriptors.getALogP();
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {

        CurOutput.setMainResultValue(Descriptors[0]);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Descriptors[0])); // LogP
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setSkipExperimental(SkipExperimentalInAD);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.9, 0.75);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.5);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
        // Sets final AD index
        ADIndexADI ADI = new ADIndexADILogP(adq.getIndexADI(), CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class), 
                CurOutput.getADIndex(ADIndexMaxError.class));
        
        CurOutput.setADI(ADI);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        double LogP_threshold_green = 3;
        double LogP_threshold_red = 8;
    
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val < LogP_threshold_green)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val < LogP_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
    }
}
