package insilico.skin_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.skin_irfmn.desciptors.EmbeddedDescriptors;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismSkinIRFMN extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismSkinIRFMN.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_skin_irfmn.xml";
        
    
    public ismSkinIRFMN() 
            throws InitFailureException {
        super(ModelData);
        
        
        // Define no. of descriptors
        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nDB";
        this.DescriptorsNames[1] = "IC1";
        this.DescriptorsNames[2] = "SIC2";
        this.DescriptorsNames[3] = "GATS8s";
        this.DescriptorsNames[4] = "nRCHO";
        this.DescriptorsNames[5] = "CATS2D_1_NL";
        this.DescriptorsNames[6] = "F4(C..O)";
        this.DescriptorsNames[7] = "MLOGP";
        
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted skin sensitization activity";
        
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

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            // nDB is fixed manually to be consistent with D7 implementation
            // NO2 groups are normalized to N(=O)=O so additional double
            // bonds are not counted properly
            double nDB = embeddedDescriptors.getNDblBo();
            int nNO2 = 0;
            if (embeddedDescriptors.getNRNO2() != Descriptor.MISSING_VALUE)
                nNO2 += embeddedDescriptors.getNRNO2();
            if (embeddedDescriptors.getNArNO2() != Descriptor.MISSING_VALUE)
                nNO2 += embeddedDescriptors.getNArNO2();
            if (nDB != Descriptor.MISSING_VALUE)
                nDB += nNO2;
            
            Descriptors = new double[DescriptorsSize];
            Descriptors[0] = nDB;
            Descriptors[1] = embeddedDescriptors.getIC1();
            Descriptors[2] = embeddedDescriptors.getSIC2();
            Descriptors[3] = embeddedDescriptors.getGATS8s();
            Descriptors[4] = embeddedDescriptors.getNRCHO();
            Descriptors[5] = embeddedDescriptors.getCATS2D_1_NL();
            Descriptors[6] = embeddedDescriptors.getF4_C_O();
            Descriptors[7] = embeddedDescriptors.getMLOGP();
                        
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        boolean CartPrediction = SkinCartTree.Predict(Descriptors[0], Descriptors[1], 
                Descriptors[2], Descriptors[3], Descriptors[4], Descriptors[5], 
                Descriptors[6], Descriptors[7]);
        
        int MainResult = (CartPrediction ? 1 : 0);
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // skin classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for skin value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        
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
