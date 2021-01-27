package insilico.rba_cerapp;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAEstrogenBindCerapp;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;


import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class ismEstrogenBindingCerapp extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_estrogen_cerapp.xml";
    
    private final String SarpyFullAlertSet;
    
    public ismEstrogenBindingCerapp() 
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        SarpyFullAlertSet = AlertEncoding.MergeAlertIds((new SAEstrogenBindCerapp()).getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted ER-mediated effect";
        this.ResultsName[1] = "No. alerts for activity";
        this.ResultsName[2] = "No. alerts for possible activity";
        this.ResultsName[3] = "No. alerts for non-activity";
        this.ResultsName[4] = "No. alerts for possible non-activity";
        
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
        blocks.add(InsilicoConstants.SA_BLOCK_ESTROGEN_BIND_CERAPP);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(SarpyFullAlertSet, a.getId())) {
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
        
        int nAct=0, nPosAct=0, nInact=0, nPosInact=0;

        try {
            AlertList FoundSAs = GetCalculatedAlert();
            for (Alert CurAlert : FoundSAs.getSAList()) {
                if (CurAlert.getBoolProperty(InsilicoConstants.KEY_ALERT_ER_ACTIVE))
                    nAct++;
                if (CurAlert.getBoolProperty(InsilicoConstants.KEY_ALERT_ER_ACTIVE_POSSIBLE))
                    nPosAct++;
                if (CurAlert.getBoolProperty(InsilicoConstants.KEY_ALERT_ER_INACTIVE))
                    nInact++;
                if (CurAlert.getBoolProperty(InsilicoConstants.KEY_ALERT_ER_INACTIVE_POSSIBLE))
                    nPosInact++;
            }
        } catch (CloneNotSupportedException ex) {
            // unable to get SAs, return error
            return MODEL_ERROR;
        }
        
        if (nAct > 0)
            MainResult = 2;
        else if (nPosAct > 0)
            MainResult = 3;
        else if (nInact > 0)
            MainResult = 0;
        else if (nPosInact > 0)
            MainResult = 1;
        else
            MainResult = -1;

        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // ER classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for ER value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Integer.toString(nAct);
        Res[2] = Integer.toString(nPosAct);
        Res[3] = Integer.toString(nInact);
        Res[4] = Integer.toString(nPosInact);
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(2);
        adq.AddMappingToNegativeValue(3);
        adq.AddMappingToNegativeValue(0);
        adq.AddMappingToNegativeValue(1);
        adq.setMoleculesForIndexSize(3);
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
        if (Val == -1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else if (Val == 2)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val == 3)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        
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
