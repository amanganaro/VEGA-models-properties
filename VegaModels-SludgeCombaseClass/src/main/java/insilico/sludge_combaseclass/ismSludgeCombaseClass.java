package insilico.sludge_combaseclass;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import protoqsar.filters.BiocideFilter;
import protoqsar.sludge.SludgeDescriptors;
import protoqsar.sludge.SludgeQualitativeTree;

import java.util.ArrayList;

/**
 *
 * @author User
 */

public class ismSludgeCombaseClass extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismSludgeCombaseClass.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_sludge_combaseClass.xml";
    
    
    public ismSludgeCombaseClass() 
            throws InitFailureException {
        super(ModelData);
        
        
        // Define no. of descriptors
        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ATSC4p";
        this.DescriptorsNames[1] = "ETA_BetaP_ns_d";
        this.DescriptorsNames[2] = "GATS1i";
        this.DescriptorsNames[3] = "GATS3c";
        this.DescriptorsNames[4] = "maxHother";
        this.DescriptorsNames[5] = "minsCH3";
        this.DescriptorsNames[6] = "minwHBa";
        this.DescriptorsNames[7] = "SpMax1_Bhm";
                
        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted toxicity activity";
        this.ResultsName[1] = "Biocide filter";
        
        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
                
    }
    
    
    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        return blocks;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            SludgeDescriptors SludgeDescEngine = new SludgeDescriptors();
            SludgeDescEngine.Calculate(this.CurMolecule);
            double[] algaeDesc = SludgeDescEngine.GetDescriptorsForSludgeClassificationModel();
            
            Descriptors = new double[DescriptorsSize];
            for (int i=0; i<DescriptorsSize; i++)
                Descriptors[i] = algaeDesc[i];
                
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        SludgeQualitativeTree SludgeTree = new SludgeQualitativeTree();
        int MainResult = SludgeTree.Predict(Descriptors);
        
        CurOutput.setMainResultValue(MainResult);
        
        String biocide = BiocideFilter.ApplyFilter(CurMolecule);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // toxicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for sludge toxicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = biocide;
        
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
        double ADIValue = adq.getIndexADI() * acfContribution;

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
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
