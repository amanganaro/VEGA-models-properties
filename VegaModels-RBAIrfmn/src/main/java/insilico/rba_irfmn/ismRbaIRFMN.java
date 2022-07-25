package insilico.rba_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;

import insilico.core.tools.utils.ModelUtilities;
import insilico.rba_irfmn.descriptors.EmbeddedDescriptors;
import lombok.extern.log4j.Log4j;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class ismRbaIRFMN extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_rba_irfmn.xml";

    private int FinalNode;
            
     
    public ismRbaIRFMN() 
            throws InitFailureException {
        super(ModelData);
        
        // Define no. of descriptors
        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "X2v";
        this.DescriptorsNames[1] = "MATS6m";
        this.DescriptorsNames[2] = "MATS8v";
        this.DescriptorsNames[3] = "MATS5p";
        this.DescriptorsNames[4] = "BEH2e";
        this.DescriptorsNames[5] = "BEH1p";
        this.DescriptorsNames[6] = "nArOH";
        this.DescriptorsNames[7] = "MLogP";
        
        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted activity";
        this.ResultsName[1] = "Classification tree final node";
        
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
            
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.getX2v();
            Descriptors[1] = embeddedDescriptors.getMATS6m();
            Descriptors[2] = embeddedDescriptors.getMATS8v();
            Descriptors[3] = embeddedDescriptors.getMATS5p();
            Descriptors[4] = embeddedDescriptors.getBEH2e();
            Descriptors[5] = embeddedDescriptors.getBEH1p();
            Descriptors[6] = embeddedDescriptors.getNArOH();
            Descriptors[7] = embeddedDescriptors.getMLogP();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        int MainResult;

        double V64 = Descriptors[0];
        double V97 = Descriptors[1];
        double V107 = Descriptors[2];
        double V119 = Descriptors[3];
        double V180 = Descriptors[4];
        double V182 = Descriptors[5];
        double V221 = Descriptors[6];
        double V250 = Descriptors[7];

        if (V221 <= 0.5) {
            if (V64 <= 4.55) {
                MainResult = 0;
                FinalNode = 4;
            } else {
                if (V250 <= 3.284) {
                    MainResult = 0;
                    FinalNode = 8;
                } else {
                    if (V119 <= -0.137) {
                        MainResult = 0;
                        FinalNode = 12;
                    } else {
                        if (V97 <= -0.011) {
                            if (V64 <= 9.887) {
                                MainResult = 1;
                                FinalNode = 18;
                            } else {
                                MainResult = 0;
                                FinalNode = 19;
                            }
                        } else {
                            if (V182 <= 3.9305) {
                                MainResult = 0;
                                FinalNode = 20;
                            } else {
                                MainResult = 1;
                                FinalNode = 21;
                            }
                        }
                    }
                }
            }
        } else {
            if (V107 <= -0.0375) {
                MainResult = 1;
                FinalNode = 6;
            } else {
                if (V180 <= 3.8035) {
                    if (V97 <= -0.4495) {
                        MainResult = 1;
                        FinalNode = 14;
                    } else {
                        MainResult = 0;
                        FinalNode = 15;
                    }
                } else {
                    MainResult = 1;
                    FinalNode = 11;
                }                                
            }
        }
    
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // RBA classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for toxicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Integer.valueOf(FinalNode).toString();
        
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
        ADI.SetThresholds(0.85, 0.7);
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
        
    }
    
}
