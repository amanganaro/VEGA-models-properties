package insilico.ld50quail;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertBlock;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.ld50quail.alerts.SALD50Quail;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismLD50Quail extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismLD50Quail.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_ld50_quail.xml";

    private final String SarpyFullAlertSet;
    private AlertList LD50FoundSAs = null;

    public ismLD50Quail()
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        SarpyFullAlertSet = AlertEncoding.MergeAlertIds((new SALD50Quail()).getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LD50 class";
        
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
        blocks.add(InsilicoConstants.SA_BLOCK_MICRONUCLEUS_INVITRO_MODEL);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        if (LD50FoundSAs == null) {
            try {
                AlertBlock ab = new SALD50Quail();
                LD50FoundSAs = ab.Calculate(CurMolecule);
            } catch (Exception e) {
                throw new CloneNotSupportedException("Unable to run alerts");
            }
        }
        return LD50FoundSAs;
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

        LD50FoundSAs = null;

        int MainResult = -1;
        
        // Check just the first alert foud (it will be the more accurate,
        // alerts are checked following their accuracy order)
        try {
            AlertList FoundSAs = GetCalculatedAlert();
            if (FoundSAs.getSAList().isEmpty())
                MainResult = -1;
            else {

                int ResL = -1;
                int ResH = -1;

                for (Alert a : FoundSAs.getSAList()) {
                    if (a.getBoolProperty("LD50QUAIL_L")) {
                        if (ResL != -1) continue;
                        ResL = a.getBoolProperty("LD50QUAIL_ACTIVE") ? 1 : 0;
                    }
                    if (a.getBoolProperty("LD50QUAIL_H")) {
                        if (ResH != -1) continue;
                        ResH = a.getBoolProperty("LD50QUAIL_ACTIVE") ? 1 : 0;
                    }
                }

                if ( (ResL==0) && (ResH==0) ) MainResult = 0;
                if ( (ResL==0) && (ResH==1) ) MainResult = 2;
                if ( (ResL==1) && (ResH==1) ) MainResult = 2;
                if ( (ResL==1) && (ResH==0) ) MainResult = 1;
                if ( (ResL==-1) && (ResH==1) ) MainResult = 2;
                if ( (ResL==0) && (ResH==-1) ) MainResult = 0;
                if ( (ResL==-1) && (ResH==0) ) MainResult = 0;
                if ( (ResL==1) && (ResH==-1) ) MainResult = 2;
                if ( (ResL==-1) && (ResH==-1) ) MainResult = -1;

            }
           
        } catch (CloneNotSupportedException ex) {
            // unable to get SAs, return error
            return MODEL_ERROR;
        }
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); //  classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for LD50 value " + MainResult);
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
        adq.AddMappingToPositiveValue(2);
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
