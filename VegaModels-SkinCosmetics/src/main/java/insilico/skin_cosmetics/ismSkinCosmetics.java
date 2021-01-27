package insilico.skin_cosmetics;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SASkinNcstox;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.skin_caesar.ismSkinCaesar;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class ismSkinCosmetics extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_skin_cosmetics.xml";
    
    private final static boolean CALCULATE_AD = false;
    
    private final String SarpyFullAlertSet;
    private ToxTreeSkinClassification ToxTreeSkin;
    private ismSkinCaesar CaesarModel;
    
    public ismSkinCosmetics() 
            throws InitFailureException {
        super(ModelData);
        
        // Set SA list
        SarpyFullAlertSet = AlertEncoding.MergeAlertIds((new SASkinNcstox()).getAlerts());
        
        // Init needed sub-models
        ToxTreeSkin = new ToxTreeSkinClassification();
        CaesarModel = new ismSkinCaesar();
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted skin sensitization activity";
        this.ResultsName[1] = "Predicted consensus reliability";
        this.ResultsName[2] = "Prediction from VEGA Caesar model";
        this.ResultsName[3] = "Prediction with ToxTree classes";
        this.ResultsName[4] = "Prediction with fragment-based model";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        
    }
    
    
    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
        blocks.add(InsilicoConstants.SA_BLOCK_SKIN_SENS_NCSTOX);
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
        
        // Calculate CAESAR and ToxTree
        InsilicoModelOutput oCaesar;
        ToxTreeSkinClassification.SkinClass oToxTree;
        try {
            oCaesar = CaesarModel.Execute(CurMolecule);
            oToxTree = ToxTreeSkin.Calculate(CurMolecule);
        } catch (GenericFailureException e) {
            return MODEL_ERROR;
        }
        
        int predToxTree = 0;
        String predToxTreeText = "NON-Sensitizer";
        if (oToxTree.getId() != ToxTreeSkinClassification.SKIN_NO_CLASS) {
            predToxTree = 1;
            predToxTreeText = "Sensitizer (" + oToxTree.getName() + ")";
        }
        
        int predCaesar = -1;
        boolean CaesarExp = false;
        String predCaesarText = "";
        if (oCaesar.getStatus() == -1) {
            predCaesarText = "N.A.";
        } else if (oCaesar.HasExperimental()) {
            predCaesarText = oCaesar.getAssessment();
            predCaesar = (int)oCaesar.getExperimental();
            CaesarExp = true;
        } else {
            predCaesarText = oCaesar.getAssessment();
            predCaesar = (int)oCaesar.getMainResultValue();            
        }
        
        
        // Check alerts
        int predSarpy = -1;
        String predSarpyText = "N.A.";
        try {            
            for (Alert a : CurMolecule.GetAlerts().getSAList())
                if (AlertEncoding.ContainsAlert(SarpyFullAlertSet, a.getId())) {
                    if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_SKIN_SENS)) {
                        predSarpy = 1;
                        predSarpyText = "Sensitizer";
                        break;
                    }
                    if (a.getBoolProperty(InsilicoConstants.KEY_ALERT_SKIN_NON_SENS)) {
                        predSarpy = 0;
                        predSarpyText = "NON-Sensitizer";
                        break;
                    }
                }                            
        } catch (Exception e) {
            return MODEL_ERROR;
        }        
        
        
        // Model - consensus
        int FinalPred = -1;
        int FinalRel = -1; // 0:exp 1:good 2:moderate 3:low
        
        // Check exp
        if (CaesarExp) {
            FinalPred = predCaesar;
            FinalRel = 0;
        } else if (predCaesar != -1) {
            
            if (oCaesar.getADI().GetAssessmentClass() == ADIndex.INDEX_HIGH) {
                
                // Caesar good rel
                
                if (predSarpy == -1) {
                    
                    // sarpy not available
                    if (predCaesar == predToxTree) {
                        FinalPred = predCaesar;
                        FinalRel = 1;
                    } else {
                        FinalPred = predCaesar;
                        FinalRel = 2;
                    }
                } else {
                    
                    // with all 3 models
                    if (predCaesar == predSarpy) {
                        FinalPred = predCaesar;
                        FinalRel = 1;
                    } else {
                        if (predSarpy == predToxTree)
                            FinalPred = predSarpy;
                        else 
                            FinalPred = predCaesar;
                        FinalRel = 2;
                    }
                    
                }
                
            } else if (oCaesar.getADI().GetAssessmentClass() == ADIndex.INDEX_MEDIUM) {
            
                // Caesar moderate rel
                double CaesarConc = oCaesar.getADIndex(ADIndexConcordance.class).GetIndexValue();
                
                if (predSarpy == -1) {
                    
                    // sarpy not available
                    FinalPred = predCaesar;
                    if (predCaesar == predToxTree) {
                        if (CaesarConc > 0.8)
                            FinalRel = 1;
                        else
                            FinalRel = 2;                            
                    } else {
                        if (CaesarConc > 0.8)
                            FinalRel = 2;
                        else
                            FinalRel = 3;                                                    
                    }
                    
                } else {
                    
                    // with all 3 models
                    if ( (predCaesar == predSarpy) && (predCaesar == predToxTree) ) {
                        
                        FinalPred = predCaesar;
                        FinalRel = 1;
                    
                    } else if (predToxTree == predSarpy) {
                        
                        FinalPred = predSarpy;
                        if ( (predCaesar == 0) && (CaesarConc > 0.8))
                            FinalRel = 3;
                        else
                            FinalRel = 2;
                        
                    } else {
                        
                        FinalPred = predCaesar;
                        if (predCaesar == 0) {
                            
                            if ( (predSarpy == 0) && (CaesarConc > 0.8) ) FinalRel = 1;
                            if ( (predSarpy == 1) && (CaesarConc > 0.8) ) FinalRel = 2;
                            if ( (predSarpy == 0) && (CaesarConc <= 0.8) ) FinalRel = 2;
                            if ( (predSarpy == 1) && (CaesarConc <= 0.8) ) FinalRel = 3;
                            
                        } else {
                            
                            FinalRel = 2;
                            
                        }
                                              
                    }
                }
                
            } else {
                
                // Caesar low rel
                if (predSarpy == -1) {
                    
                    // sarpy not available
                    if (predCaesar == predToxTree) {
                        FinalPred = predCaesar;
                        FinalRel = 2;
                    } else {
                        FinalPred = predToxTree;
                        FinalRel = 3;
                    }
                } else {
                    
                    // with all 3 models
                    if ( (predCaesar == predSarpy) && (predCaesar == predToxTree) ) {
                        FinalPred = predCaesar;
                        FinalRel = 1;
                    } else {
                        if (predSarpy == predToxTree)
                            FinalPred = predSarpy;
                        else 
                            FinalPred = predCaesar;
                        FinalRel = 2;
                    }
                    
                }
                
                
            }
        } 
        
        String FinalRelTxt = "";
        switch (FinalRel) {
            case 0: FinalRelTxt = "EXPERIMENTAL"; break;
            case 1: FinalRelTxt = "Good reliability"; break;
            case 2: FinalRelTxt = "Moderate reliability"; break;
            case 3: FinalRelTxt = "Low reliability"; break;
            default: FinalRelTxt = "N.A."; 
        }
        
        int MainResult = FinalPred;
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // skin classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for skin value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = FinalRelTxt;
        Res[2] = predCaesarText;
        Res[3] = predToxTreeText;
        Res[4] = predSarpyText;
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // for now - AD not calculated
        if (!CALCULATE_AD)
            return InsilicoModel.AD_ERROR;
        
        
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
