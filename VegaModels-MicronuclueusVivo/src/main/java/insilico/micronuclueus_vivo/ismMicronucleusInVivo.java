package insilico.micronuclueus_vivo;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAMicronucleusInVivo;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.knn.insilicoKnnQualitative;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;


import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class ismMicronucleusInVivo extends InsilicoModel {
    
    private static final short PRED_NA = 0;
    private static final short PRED_AGREEMENT = 1;
    private static final short PRED_NOT_AGREEMENT = -1;
    private static final short PRED_ONLY_ON_ONE = -2;
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_micronucleus_vivo.xml";
    private static final String TS_KNN_PATH = "/data/ts_micronucleus_vivo_knn.dat";
    
    private final String SarpyFullAlertSet;
    private TrainingSet TS_KNN;
    private insilicoKnnQualitative KNN;
    private insilicoKnnPrediction KnnPrediction;

    public boolean KnnSkipExperimental;
    private boolean KNN_TS_ISBUILDING;

    private void BuildKNNAndModel() {

        // Build KNN model object
        KNN =  new insilicoKnnQualitative();
        KNN.setNeighboursNumber(4);
        KNN.setMinSimilarity(0.75);
        KNN.setMinSimilarityForSingleResult(0.85);
        KNN.setEnhanceWeightFactor(1);
        KNN.setUseExperimentalRange(true);

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 7;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted activity";
        this.ResultsName[1] = "Consensus status";
        this.ResultsName[2] = "Sarpy model result";
        this.ResultsName[3] = "Sarpy used alerts reliability";
        this.ResultsName[4] = "KNN model result";
        this.ResultsName[5] = "KNN reliability";
        this.ResultsName[6] = "KNN molecules used for prediction";


        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

    }
    
    public ismMicronucleusInVivo()
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        SarpyFullAlertSet = AlertEncoding.MergeAlertIds((new SAMicronucleusInVivo()).getAlerts());


        // Read TS for the KNN model
        try {
            URL u = getClass().getResource(TS_KNN_PATH);
            ObjectInputStream in = new ObjectInputStream(u.openStream());
            TS_KNN = (TrainingSet) in.readObject();
            in.close();
            this.KNN_TS_ISBUILDING = false;
        } catch (IOException | ClassNotFoundException | NullPointerException e) {
            this.ProcessKNNTrainingSet();
            log.info("KNN TRAINING SET PROCESSED");
//            throw new InitFailureException("Unable to load training set for KNN model - " + e.getMessage());
        }

//
        KNN =  new insilicoKnnQualitative();
        KNN.setNeighboursNumber(4);
        KNN.setMinSimilarity(0.75);
        KNN.setMinSimilarityForSingleResult(0.85);
        KNN.setEnhanceWeightFactor(1);
        KNN.setUseExperimentalRange(false);
        KnnSkipExperimental = false;

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 7;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted activity";
        this.ResultsName[1] = "Consensus status";
        this.ResultsName[2] = "Sarpy model result";
        this.ResultsName[3] = "Sarpy used alerts reliability";
        this.ResultsName[4] = "KNN model result";
        this.ResultsName[5] = "KNN reliability";
        this.ResultsName[6] = "KNN molecules used for prediction";


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
        blocks.add(InsilicoConstants.SA_BLOCK_MICRONUCLEUS_INVIVO);
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
        short FinalPredStatus = PRED_NA;
        
        //// Sub-model KNN
        
        int KNNResult = -1;
        int KNNReliability = -1;

        if(KNN_TS_ISBUILDING) {
            BuildKNNAndModel();
        }

        try {
            KnnPrediction = KNN.Calculate(CurMolecule, TS_KNN, KnnSkipExperimental);
            if (KnnPrediction.getStatus() <= insilicoKnnPrediction.KNN_ERROR)
                KNNResult = -1;
            else
                KNNResult = (int)KnnPrediction.getPrediction();
            
            InsilicoModelOutput KnnOut = new InsilicoModelOutput();
            KnnOut.setMainResultValue(KNNResult);
            ADIndexADI KnnADI = CalculateADForKNN(KnnOut);
            if (KnnADI != null) {
                KNNReliability = KnnADI.GetAssessmentClass();
            }
            
        } catch (GenericFailureException ex) {
            return MODEL_ERROR;
        }
                
        
        //// Sub-model Sarpy
        
        int step_1_nActive = 0;
        int step_1_nInactive = 0;
        int step_2_nActive = 0;
        int step_2_nInactive = 0;
        int step_3_nAlerts = 0;

        try {            
            for (Alert a : CurMolecule.GetAlerts().getSAList())
                if (AlertEncoding.ContainsAlert(SarpyFullAlertSet, a.getId())) {

                    if (a.getNumericProperty(InsilicoConstants.KEY_ALERT_MICRONUCLEUS_INVIVO_SA_BLOCK) == SAMicronucleusInVivo.SA_REL_HIGH) {
                        if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_MICRONUCLEUS_ACTIVE))
                            step_1_nActive++;
                        else
                            step_1_nInactive++;
                    }

                    if (a.getNumericProperty(InsilicoConstants.KEY_ALERT_MICRONUCLEUS_INVIVO_SA_BLOCK) == SAMicronucleusInVivo.SA_REL_MODERATE) {
                        if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_MICRONUCLEUS_ACTIVE))
                            step_2_nActive++;
                        else
                            step_2_nInactive++;
                    }

                    if (a.getNumericProperty(InsilicoConstants.KEY_ALERT_MICRONUCLEUS_INVIVO_SA_BLOCK) == SAMicronucleusInVivo.SA_REL_LOW) {
                        step_3_nAlerts++; // all SA in block 3 are active
                    }
                }            
        } catch (Exception e) {
            return MODEL_ERROR;
        }
        
        int SarpyResult = -1;
        int SarpyReliability = -1;
        
        if (step_1_nActive + step_1_nInactive > 0) {
            
            // Step 1 (block with high reliable SAs)
            
            SarpyReliability = 3;
            if (step_1_nActive >= step_1_nInactive)
                SarpyResult = 1;
            else
                SarpyResult = 0;
            
        } else {
            
            if (step_2_nActive + step_2_nInactive > 0) {
            
                // Step 2 (block with moderate reliable SAs)
            
                SarpyReliability = 2;
                if (step_2_nActive > 0)
                    SarpyResult = 1;
                else
                    SarpyResult = 0;   

            } else {
    
                if (step_3_nAlerts > 0) {
                
                    // Step 3 (block with low reliable SAs, all active SAs)
            
                    SarpyReliability = 1;
                    SarpyResult = 1;
                    
                }                
            }
        }
        
        
        // Consensus
        
        if (KNNResult == -1) {
            
            FinalPredStatus = PRED_ONLY_ON_ONE;
            switch (SarpyResult) {
                case 0:
                    MainResult = 0; break;
                case 1: 
                    MainResult = 1; break;
                default: {
                    FinalPredStatus = PRED_NA;
                    MainResult = -1;
                }
            }
            
        } else if (SarpyResult == -1) {
            
            FinalPredStatus = PRED_ONLY_ON_ONE;
            switch (KNNResult) {
                case 0:
                    MainResult = 0; break;
                case 1: 
                    MainResult = 1; break;
                default: {
                    FinalPredStatus = PRED_NA;
                    MainResult = -1;
                }
            }
            
        } else {
            
            if ( (SarpyResult == 1) && (KNNResult == 1) ) {
                FinalPredStatus = PRED_AGREEMENT;
                MainResult = 1;
            } else if ( (SarpyResult == 0) && (KNNResult == 0) ) {
                FinalPredStatus = PRED_AGREEMENT;
                MainResult = 0;
            } else {
                // no agreement
                FinalPredStatus = PRED_NOT_AGREEMENT;
                if (SarpyReliability == KNNReliability) {
                    
                    MainResult = -1;
                
                } else if (SarpyReliability > KNNReliability) {
                    
                    if (SarpyReliability > 1) {
                        if (SarpyResult == 1)
                            MainResult = 1;
                        else if (SarpyResult == 0)
                            MainResult = 0;
                    } else
                        MainResult = -1;
                    
                } else if (KNNReliability > SarpyReliability) {
                    
                    if (KNNReliability > 1) {
                        if (KNNResult == 1)
                            MainResult = 1;
                        else if (KNNResult == 0)
                            MainResult = 0;
                    } else
                        MainResult = -1;

                }

            }
        }
        

        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // mn classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for mn value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        
        switch (FinalPredStatus) {
            case PRED_AGREEMENT:
                Res[1] = "Both predictions in agreement"; break;
            case PRED_NOT_AGREEMENT:
                Res[1] = "Predictions NOT in agreement"; break;
            case PRED_ONLY_ON_ONE:
                Res[1] = "Based only on ONE available prediction"; break;
            default:
                Res[1] = "Not available";
        }
        
        try {
            Res[2] = this.GetTrainingSet().getClassLabel(SarpyResult); // sarpy classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for mn value " + MainResult);
            Res[2] = Integer.toString(SarpyResult);
        }
        String SarpyRel = "-";
        if (SarpyReliability==3) SarpyRel = "High";
        if (SarpyReliability==2) SarpyRel = "Moderate";
        if (SarpyReliability==1) SarpyRel = "Low";
        Res[3] = SarpyRel;
        try {
            Res[4] = this.GetTrainingSet().getClassLabel(KNNResult); // mn classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for mn value " + MainResult);
            Res[4] = Integer.toString(KNNResult);
        }
        String KnnRel = "-";
        if (KNNReliability== ADIndex.INDEX_HIGH) KnnRel = "High";
        if (KNNReliability== ADIndex.INDEX_MEDIUM) KnnRel = "Moderate";
        if (KNNReliability== ADIndex.INDEX_LOW) KnnRel = "Low";
        Res[5] = String.valueOf(KnnRel);
        Res[6] = String.valueOf(KnnPrediction.getNeighbours().size());
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    //todo
    private ADIndexADI CalculateADForKNN(InsilicoModelOutput KnnOutput) {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS_KNN);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        // Ad is performed on the K molecules used for the KNN model
        adq.setMoleculesForIndexSize(KnnPrediction.getNeighbours().size());

        // Nothing calculated if KNN gave no prediction
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == Descriptor.MISSING_VALUE)
            return null;
        
        if (!adq.Calculate(CurMolecule, KnnOutput))
            return null;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)KnnOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.65);
            ((ADIndexAccuracy)KnnOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.8, 0.65);
            ((ADIndexConcordance)KnnOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.8, 0.65);
        } catch (Throwable e) {
            return null;
        }
        
        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS_KNN);
        if (!adacf.Calculate(CurMolecule, KnnOutput))
            return null;
        
        // Sets final AD index
        double acfContribution = KnnOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.8, 0.65);
        KnnOutput.setADI(ADI);
        
        return ADI;        
    }
    
    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.AddMappingToNegativeValue(-1);
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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.9, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.9, 0.5);
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
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
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
    
    // todo
    public void ProcessKNNTrainingSet() throws InitFailureException {
        try {
            this.setSkipADandTSLoading(false);
            this.KnnSkipExperimental = true;
            this.KNN_TS_ISBUILDING = true;
            TS_KNN = new TrainingSet();
            String TSPath = "/data/ts_micronucleus_vivo_knn.dat";
            String[] buf = TSPath.split("/");
            String DatName = buf[buf.length-1];
            TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
            TS_KNN.Build(TSPath, this, true, false);
            TS_KNN.SerializeToFile("out_ts/" + DatName);
            File sourceFile = new File("out_ts/" + DatName);
            File destinationFile = new File("VegaModels-MicronuclueusVivo\\src\\main\\resources\\data\\ts_micronucleus_vivo_knn.dat");
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
            throw  new InitFailureException((ex));
        }
    }
    
//    @Override
//    public void ProcessTrainingSet() throws Exception {
//        this.setSkipADandTSLoading(false);
//
//        // When building TS, KNN prediction is run skipping exp (LOO)
//        // so to have unbiased predictions
//        System.out.println("Skipping exp values in KNN model");
//        this.KnnSkipExperimental = true;
//
//        TrainingSet TS = new TrainingSet();
//        String TSPath = this.getInfo().getTrainingSetURL();
//        String[] buf = TSPath.split("/");
//        String DatName = buf[buf.length-1];
//        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
//        TS.Build(TSPath, this);
//        TS.SerializeToFile(DatName);
//    }


    
    
    
}
