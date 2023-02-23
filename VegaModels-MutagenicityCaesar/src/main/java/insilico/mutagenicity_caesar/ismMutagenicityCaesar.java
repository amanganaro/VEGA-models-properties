package insilico.mutagenicity_caesar;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.mutagenicity_caesar.descriptors.EmbeddedDescriptors;
import insilico.mutagenicity_caesar.descriptors.MutagenDescriptors;
import insilico.mutagenicity_caesar.descriptors.weights.EState;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismMutagenicityCaesar extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismMutagenicityCaesar.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_muta_caesar.xml";
    
    private final ModelSVM SVM;
    private final String[] MutagenSA;
    private final String[] SuspectMutagenSA;
    private final String MutagenModelSA;

    
    public ismMutagenicityCaesar() 
            throws InitFailureException {
        super(ModelData);
        
        // Builds model objects
        URL uModel, uRanges;
        uModel = getClass().getResource("/data/muta25_SVM8_16.model");
        uRanges = getClass().getResource("/data/muta_ranges.csv");
        try {
            SVM = new ModelSVM(uModel.openStream(), uRanges.openStream());
        } catch (IOException ex) {
            throw new InitFailureException("Unable to initialize SVM");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 25;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "SsCH3";
        this.DescriptorsNames[1] = "SdCH2";
        this.DescriptorsNames[2] = "SssCH2";
        this.DescriptorsNames[3] = "SdsCH";
        this.DescriptorsNames[4] = "SaaCH";
        this.DescriptorsNames[5] = "SsssCH";
        this.DescriptorsNames[6] = "SdssC";
        this.DescriptorsNames[7] = "SaasC";
        this.DescriptorsNames[8] = "SaaaC";
        this.DescriptorsNames[9] = "SssssC";
        this.DescriptorsNames[10] = "SsNH2";
        this.DescriptorsNames[11] = "StN";
        this.DescriptorsNames[12] = "SdsN";
        this.DescriptorsNames[13] = "SaaN";
        this.DescriptorsNames[14] = "SsssN";
        this.DescriptorsNames[15] = "SsaaN";
        this.DescriptorsNames[16] = "SsOH";
        this.DescriptorsNames[17] = "SdO";
        this.DescriptorsNames[18] = "SssO";
        this.DescriptorsNames[19] = "SaaO";
        this.DescriptorsNames[20] = "SHCHnX";
        this.DescriptorsNames[21] = "Gmin";
        this.DescriptorsNames[22] = "idwbar";
        this.DescriptorsNames[23] = "nrings";
        this.DescriptorsNames[24] = "ALogP";
        
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Mutagen activity";
        
        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        
        // Builds alerts list for second step of the model:
        // SA1, SA6, SA12, SA13, SA14, SA16, SA18, SA21,
        // SA22, SA25, SA28b, SA29, 
        int[] MutagenAlertsId = {0,5,11,12,13,15,17,20,21,24,28,30};
        MutagenSA = new String[MutagenAlertsId.length];
        for (int i=0; i<MutagenAlertsId.length; i++)
            MutagenSA[i] = AlertEncoding.BuildAlertId(InsilicoConstants.SA_BLOCK_MUTAGEN_BENIGNI_BOSSA, (MutagenAlertsId[i]+1));

        // Builds alerts list for third step of the model:
        // SA7, SA8, SA19, SA27
        int[] SuspectMutagenAlertsId = {6,7,18,26};
        SuspectMutagenSA = new String[SuspectMutagenAlertsId.length];
        for (int i=0; i<SuspectMutagenAlertsId.length; i++)
            SuspectMutagenSA[i] = AlertEncoding.BuildAlertId(InsilicoConstants.SA_BLOCK_MUTAGEN_BENIGNI_BOSSA, (SuspectMutagenAlertsId[i]+1));

        // Build list of all alerts used by the model
        String AllSA = "";
        for (String s : MutagenSA)
            AllSA = AlertEncoding.MergeAlertIds(AllSA, s);
        for (String s : SuspectMutagenSA)
            AllSA = AlertEncoding.MergeAlertIds(AllSA, s);
        this.MutagenModelSA = AllSA;
    }
    

    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
        blocks.add(InsilicoConstants.SA_BLOCK_MUTAGEN_BENIGNI_BOSSA);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(MutagenModelSA, a.getId())) {
                FoundSAs.add((Alert)a.clone());
            }
        return FoundSAs;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
            
            // Calculates original open Mutagenicity descriptors
            double[] MutaDescriptors = MutagenDescriptors.Calculate(CurMolecule);

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            // Fragment counter descriptors (0-20)
            for (int i=0; i<21; i++)
                Descriptors[i] = MutaDescriptors[i];

            Descriptors[21] = embeddedDescriptors.getGmin();
            Descriptors[22] = embeddedDescriptors.getIDWBAR();
            Descriptors[23] = CurMolecule.GetSSSR().getAtomContainerCount(); // nRings calculated directly here
            Descriptors[24] = embeddedDescriptors.getALogP();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        int MainResult = -1;
        
        // First step: SVM

        boolean SVMResult = SVM.Calculate_Prediction(Descriptors);

        if (SVMResult == true) {

            MainResult = 1; // Mutagen

        } else {

            //// Second step: SA CHECK 1
            //// SA bound to mutagen activity

            String curAlerts = AlertEncoding.MergeAlertIds(CurMolecule.GetAlerts());
            
            boolean foundMuta = false;
            for (String SA : MutagenSA)
                if (AlertEncoding.ContainsAlert(curAlerts, SA)) {
                    foundMuta = true;
                    break;
                }
                
            if (foundMuta == true) {

                MainResult = 1; // Mutagen

            } else {

                //// second step: SA CHECK 2
                //// SA bound to suspect mutagen activity

                boolean foundSuspect = false;
                for (String SA : SuspectMutagenSA)
                    if (AlertEncoding.ContainsAlert(curAlerts, SA)) {
                        foundSuspect = true;
                        break;
                    }
                
                if (foundSuspect == true) {

                    MainResult = 2; // Suspect mutagen

                } else {

                    MainResult = 0; // Non mutagen

                }
            }
        }        
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // Mutagenicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for mutagenicity value " + MainResult);
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
        adq.AddMappingToPositiveValue(2); // DA VEDERE!
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(3);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.9, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.9, 0.5);
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
        ADI.SetThresholds(0.9, 0.7);
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
        else if (Val == 2)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        
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
