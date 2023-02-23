package insilico.devtox_pg;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.devtox_pg.library.DARTLibrary;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismDevToxPG extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismDevToxPG.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_devtox_pg.xml";
    
    private final DARTLibrary DevToxLib;
    private DARTLibrary.DARTResult ResCategory;
    
    
    public ismDevToxPG() 
            throws InitFailureException {
        super(ModelData);
        
        // Set model object
        DevToxLib = new DARTLibrary();
        ResCategory = null;
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Developmental and Reproductive activity";
        this.ResultsName[1] = "Predicted category";
        this.ResultsName[2] = "Predicted category description";
        this.ResultsName[3] = "Predicted sub-category";
        this.ResultsName[4] = "Matching rule/virtual compound";
        
        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();
        
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

        int MainResult = 0;

        try {
            ResCategory = DevToxLib.Check(CurMolecule);
        } catch (Throwable ex) {
            return MODEL_ERROR;
        }
        
        if (ResCategory == null) {
            MainResult = 0;
        } else {
            if (ResCategory.Index == 26) 
                MainResult = 0;
            else
                MainResult = 1;
        }
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // devtox classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for dev tox value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = ResCategory == null ? "-" : "" + ResCategory.Index;
        Res[2] = ResCategory == null ? "-" : DARTLibrary.getCategoryName(ResCategory.Index);
        Res[3] = ResCategory == null ? "-" : ResCategory.Index+ResCategory.SubIndex;
        Res[4] = ResCategory == null ? "-" : ResCategory.Structure;
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        
        // for predicted values
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        // for data in training set
        adq.AddMappingToPositiveValue(2);
        adq.AddMappingToPositiveValue(3);
        adq.AddMappingToPositiveValue(4);
        adq.AddMappingToPositiveValue(5);
        adq.AddMappingToPositiveValue(8);
        adq.AddMappingToNegativeValue(6);
        adq.AddMappingToNegativeValue(7);
        adq.AddMappingToNegativeValue(9);
        adq.AddMappingToNegativeValue(10);

        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.75);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.85, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.85, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.9, 0.65);
        CurOutput.setADI(ADI);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if ((Val == 0)||(Val == 6)||(Val == 7)||(Val == 9))
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if ((Val == 2)||(Val == 3)||(Val == 4)||(Val == 5)||(Val == 8)||(Val == 1))
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        
        // Additional assessment when structure has been found
        if (ResCategory != null) {
            String category = "\nThe molecule matches with a category from the library: " +
                    ResCategory.Index + ResCategory.SubIndex;
            CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() + category);
        }
    }
    
}
