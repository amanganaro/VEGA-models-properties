package insilico.readybio_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAReadyBioIRFMN;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
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
public class ismReadyBioIRFMN extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_readybio_irfmn.xml";
    
    private final String ReadyBioFullAlertSet;
    
    
    public ismReadyBioIRFMN() 
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        ReadyBioFullAlertSet = AlertEncoding.MergeAlertIds((new SAReadyBioIRFMN()).getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted RB activity";
        this.ResultsName[1] = "No. alerts for non RB";
        this.ResultsName[2] = "No. alerts for possible non RB";
        this.ResultsName[3] = "No. alerts for RB";
        this.ResultsName[4] = "No. alerts for possible RB";
        
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
        blocks.add(InsilicoConstants.SA_BLOCK_READY_BIO_IRFMN);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(ReadyBioFullAlertSet, a.getId())) {
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
        
        int nRB = 0, nRBpos=0, nNonRB=0, nNonRBpos=0;

        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(ReadyBioFullAlertSet, a.getId())) {
                if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_READY_BIO_NON_RB))
                    nNonRB++;
                else if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_READY_BIO_NON_RB_POSSIBLE))
                    nNonRBpos++;
                if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_READY_BIO_RB))
                    nRB++;
                else if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_READY_BIO_RB_POSSIBLE))
                    nRBpos++;
            }
        
        if (nNonRB > 0)
            MainResult = 1; // non RB
        else if (nNonRBpos > 0)
            MainResult = 2; // probable non RB
        else if (nRB > 0)
            MainResult = 0; // RB
        else if (nRBpos > 0)
            MainResult = 3; // probable RB
        else
            MainResult = -1; // unknown
        
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // RB classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for RB value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Integer.toString(nNonRB);
        Res[2] = Integer.toString(nNonRBpos);
        Res[3] = Integer.toString(nRB);
        Res[4] = Integer.toString(nRBpos);
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToPositiveValue(2);
        adq.AddMappingToNegativeValue(0);
        adq.AddMappingToNegativeValue(3);
        adq.setMoleculesForIndexSize(3);

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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.8, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.8, 0.5);
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
        else if (Val == 3)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val == 2)
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
