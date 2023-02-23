package insilico.hepatoxicty_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAHepatotoxicity;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismHepatotoxicityIrfmn extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismHepatotoxicityIrfmn.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_hepa_irfmn.xml";
    
    private final String FullAlertSet;
    
    public ismHepatotoxicityIrfmn() 
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        FullAlertSet = AlertEncoding.MergeAlertIds((new SAHepatotoxicity()).getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Hepatotoxicity";
        this.ResultsName[1] = "No. alerts for toxicity";
        this.ResultsName[2] = "No. alerts for non-toxicity";
        
        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();
        
    }
    
    
    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
        blocks.add(InsilicoConstants.SA_BLOCK_HEPATOTOXICITY);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(FullAlertSet, a.getId())) {
                FoundSAs.add((Alert)a.clone());
            }
        return FoundSAs;
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
        
        int MainResult = -1;
        
        int nTox = 0;
        int nNonTox = 0;

        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(FullAlertSet, a.getId())) {
                if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_HEPA_TOXIC))
                    nTox++;
                else if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_HEPA_NONTOXIC))
                    nNonTox++;
            }
        
        int nTot = nTox+nNonTox;
        
        if (nTot == 0)
            // no SAs
            MainResult = -1;
        else if (nTot == 1) {
            // only one SA
            if (nTox == 1)
                MainResult = 1;
            else
                MainResult = 0;
        } else {
            // more than one SA
            if (nNonTox > nTox)
                MainResult = 0;
            else
                MainResult = 1;                
        }
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // hepatotoxicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for toxicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Integer.toString(nTox);
        Res[2] = Integer.toString(nNonTox);
        
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
        
        // Sets SA
        try {
            ADCheckSA ADSA = new ADCheckSA(TS);
            if (!ADSA.Calculate(CurMolecule, GetCalculatedAlert(), CurOutput, adq.GetSimilarMolecules()))
                return InsilicoModel.AD_ERROR;
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
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
        
        // Additional assessment when some SAs are found
        try {
            AlertList curSA = GetCalculatedAlert();
            if (curSA.size() > 0) {
                String SAs = "\nThe following relevant fragments have been found: " +
                        ModelUtilities.BuildSANameList(curSA.getSAList());
                CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() + SAs);
            }
        } catch (CloneNotSupportedException e) {
            log.warn("unable to add SAs in the assessment");
        }
    }
    
}
