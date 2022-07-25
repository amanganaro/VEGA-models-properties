package insilico.carcinogenicity_sfi_classification;

import insilico.carcinogenicity_sfi_classification.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
//import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;

/**
 *
 * @author User
 */
@Log4j
public class ismCarcinogenicitySFIClassification extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_carc_sficlassification.xml";
    
    
    
    public ismCarcinogenicitySFIClassification() 
            throws InitFailureException {
        super(ModelData);
        
        // Define no. of descriptors
        this.DescriptorsSize = 9;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0]="piPC10";
        this.DescriptorsNames[1]="ATSC2p";
        this.DescriptorsNames[2]="EEig15bo";
        this.DescriptorsNames[3]="F1(N..N)";
        this.DescriptorsNames[4]="BEL3m";
        this.DescriptorsNames[5]="GATS1v";
        this.DescriptorsNames[6]="PW4";
        this.DescriptorsNames[7]="MATS3s";
        this.DescriptorsNames[8]="BEH4e";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Inhalation Carcinogenic class";
        
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

            Descriptors[0] = embeddedDescriptors.piPC10;
            Descriptors[1] = embeddedDescriptors.ATSC2p;
            Descriptors[2] = embeddedDescriptors.EEig15bo;
            Descriptors[3] = embeddedDescriptors.F1N_N;
            Descriptors[4] = embeddedDescriptors.BEL3m;
            Descriptors[5] = embeddedDescriptors.GATS1v;
            Descriptors[6] = embeddedDescriptors.PW4;
            Descriptors[7] = embeddedDescriptors.MATS3s;
            Descriptors[8] = embeddedDescriptors.BEH4e;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        int Prediction;
        
        // Tree from CART results in R
        double piPC10 = Descriptors[0];
        double ATSC2p = Descriptors[1];
        double EEig15bo = Descriptors[2];
        double F1_N_N = Descriptors[3];
        double BEL3m = Descriptors[4];
        double GATS1v = Descriptors[5];
        double PW4 = Descriptors[6];
        double MATS3s = Descriptors[7];
        double BEH4e = Descriptors[8];
        
        if (piPC10 < 6.134) {
            if (F1_N_N < 0.5) {
                if (BEL3m < 1.706) {
                    if (PW4 < 0.231) {
                        if (MATS3s < -0.0735) {
                            if (MATS3s >= -1.04) {
                                Prediction = 0;
                            } else {
                                Prediction = 1;
                            }
                        } else {
                            if (BEL3m >= -0.5935) {
                                Prediction = 0;
                            } else {
                                if (BEH4e >= 2.989) {
                                    Prediction = 0;
                                } else {
                                    Prediction = 1;
                                }
                            }
                        }
                    } else {
                        Prediction = 1;
                    }
                } else {
                    if (GATS1v >= 0.579) {
                        Prediction = 0;
                    } else {
                        Prediction = 1;
                    }
                }
            } else {
                Prediction = 1;
            }
        } else {
           if (ATSC2p >= 0.628) {
               if (piPC10 < 7.204) {
                   if (EEig15bo >= -1.612) {
                       Prediction = 0;
                   } else {
                       Prediction = 1;
                   }
               } else {
                   Prediction = 1;
               }
           } else {
               Prediction = 1;
           }
        }
        
        CurOutput.setMainResultValue(Prediction);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction); // Carc classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for carcinogenicity value " + Prediction);
            Res[0] = Integer.toString(Prediction);
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
        else 
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);

    }
}
