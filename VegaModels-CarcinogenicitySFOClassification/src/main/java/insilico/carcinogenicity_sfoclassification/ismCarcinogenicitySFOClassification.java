package insilico.carcinogenicity_sfoclassification;

import insilico.carcinogenicity_sfoclassification.descriptors.EmbeddedDescriptors;
import insilico.carcinogenicity_sfoclassification.descriptors.weights.EState;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;

/**
 *
 * @author User
 */

public class ismCarcinogenicitySFOClassification extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismCarcinogenicitySFOClassification.class);

    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_carc_sfoclassification.xml";
    
    
    
    public ismCarcinogenicitySFOClassification() 
            throws InitFailureException {
        super(ModelData);
        
        // Define no. of descriptors
        this.DescriptorsSize = 7;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0]="nS";
        this.DescriptorsNames[1]="nCIC";
        this.DescriptorsNames[2]="ATSC6s";
        this.DescriptorsNames[3]="P_VSA_logp_6";
        this.DescriptorsNames[4]="SpMaxFun";
        this.DescriptorsNames[5]="B02[C-N]";
        this.DescriptorsNames[6]="B09[C-F]";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Oral Carcinogenic class";
        
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
            
            // nCIC
            int nCIC = CurMolecule.GetSSSR().getAtomContainerCount();
            
            // Frequency and binary atom pairs
            int nSK = CurMolecule.GetStructure().getAtomCount();
            int[][] TopoMat = CurMolecule.GetMatrixTopologicalDistance();

            int B02_C_N = 0;
            int B09_C_F = 0;
            int nS = 0;

            for (int i=0; i<(nSK); i++) 
                if ((CurMolecule.GetStructure().getAtom(i).getSymbol().equalsIgnoreCase("S")))
                    nS++;
                
            for (int i=0; i<(nSK-1); i++) {
                
                for (int j=i+1; j<nSK; j++) {
                    String a = CurMolecule.GetStructure().getAtom(i).getSymbol();
                    String b = CurMolecule.GetStructure().getAtom(j).getSymbol();

                    if (TopoMat[i][j] == 2) {
                        if ( ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("N"))) || ((b.equalsIgnoreCase("C")) && (a.equalsIgnoreCase("N"))) )
                            B02_C_N = 1;
                    }
                    if (TopoMat[i][j] == 9) {
                        if ( ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("F"))) || ((b.equalsIgnoreCase("C")) && (a.equalsIgnoreCase("F"))) )
                            B09_C_F = 1;
                    }
                }
            }
            
            Descriptors[0] = nS;
            Descriptors[1] = nCIC;
            Descriptors[2] = embeddedDescriptors.ATSC6s;
            Descriptors[3] = embeddedDescriptors.P_VSA_logp_6;
            Descriptors[4] = embeddedDescriptors.SpMaxdm;
            Descriptors[5] = B02_C_N;
            Descriptors[6] = B09_C_F;            
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        int Prediction;
        
        // Tree from CART results in R
        double nS = Descriptors[0];
        double nCIC = Descriptors[1];
        double ATSC = Descriptors[2];
        double PVSA = Descriptors[3];
        double SpMax = Descriptors[4];
        double B02_C_N = Descriptors[5];
        Double B09_C_F = Descriptors[6];
        
        if (nCIC < 1.5) {
           if (PVSA < 18.59) {
               if (PVSA >= 0.4888) {
                   Prediction = 0;
               } else {
                   if (B02_C_N < 0.5) {
                       if (ATSC >= 0.02083) {
                           if (ATSC < 2.239) {
                               Prediction = 0;
                           } else {
                               if (SpMax >= 0.2) {
                                   Prediction = 0;
                               } else {
                                   if (nS >= 0.5)
                                       Prediction = 0;
                                   else
                                       Prediction = 1;
                               }
                           }
                       } else {
                           if (SpMax < (3.95 * 0.000000001)) {
                               Prediction = 0;
                           } else {
                               if (SpMax >= 1.335)
                                   Prediction = 0;
                               else
                                   Prediction = 1;
                           }
                       }
                   } else {
                       if (ATSC >= 37.71) {
                           if (ATSC >= 118.9) {
                               Prediction = 0;
                           } else {
                               if (ATSC < 102.7)
                                   Prediction = 0;
                               else
                                   Prediction = 1;
                           }
                       } else {
                           if (ATSC < 1.307) {
                               if (ATSC >= 0.8071) {
                                   Prediction = 0;
                               } else {
                                   if (SpMax >= 2.171)
                                       Prediction = 0;
                                   else
                                       Prediction = 1;
                               }
                           } else {
                               Prediction = 1;
                           }
                       }
                   }
               }
           } else {
               if (SpMax >= 3.033)
                   Prediction = 0;
               else
                   Prediction = 1;
           }
        } else {
            if (SpMax >= 1.476) {
                if (nCIC < 3.5) {
                    if (ATSC >= 68.98) {
                        Prediction = 0;
                    } else {
                        if (ATSC < 59.48) {
                            if (SpMax < 2.565) {
                                if (SpMax >= 1.918) {
                                    Prediction = 0;
                                } else {
                                    if (SpMax < 1.669)
                                        Prediction = 0;
                                    else
                                        Prediction = 1;
                                }
                            } else {
                                if (SpMax >= 2.881)
                                    Prediction = 0;
                                else
                                    Prediction = 1;
                            }
                        } else {
                            Prediction = 1;
                        }
                    }
                } else {
                    Prediction = 1;
                }
            } else {
                if (ATSC < 4.823) {
                    if (B02_C_N < 0.5)
                        Prediction = 0;
                    else
                        Prediction = 1;
                } else {
                    Prediction = 1;
                }
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
