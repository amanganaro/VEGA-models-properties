package insilico.mutagenicity_consensus;

import insilico.core.ad.item.ADIndex;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelConsensus;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.runner.InsilicoModelWrapper;
import insilico.mutagenicity_bb.ismMutagenicityBB;
import insilico.mutagenicity_caesar.ismMutagenicityCaesar;
import insilico.mutagenicity_knn.ismMutagenicityKnn;
import insilico.mutagenicity_sarpy.ismMutagenicitySarpy;


import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ismcMutagenicity extends InsilicoModelConsensus {

    private static final String ModelData = "/data/model_muta_consensus.xml";
    
    
    public ismcMutagenicity() throws InitFailureException {
        super(ModelData);
        
        // Defines results
        this.ResultsSize = 7;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Consensus Mutagen activity";
        this.ResultsName[1] = "Mutagenic Score";
        this.ResultsName[2] = "Non-Mutagenic Score";
        this.ResultsName[3] = "Model Caesar assessment";
        this.ResultsName[4] = "Model ISS assessment";
        this.ResultsName[5] = "Model SarPy assessment";
        this.ResultsName[6] = "Model KNN assessment";
    }

    
    private double GetPredictionWeight(InsilicoModelOutput Result) {
        
        if (Result.getADI().GetAssessmentClass() == ADIndex.INDEX_HIGH)
            return 0.9;
        if (Result.getADI().GetAssessmentClass() == ADIndex.INDEX_MEDIUM)
            return 0.6;
        if (Result.getADI().GetAssessmentClass() == ADIndex.INDEX_LOW)
            return 0.2;
        
        return 0;
    }
    

    @Override
    protected short CalculateModel() {

        double weightMuta = 0;
        double weightNonMuta = 0;
        int usedModels = 0;

        double weightMutaExp = 0;
        double weightNonMutaExp = 0;
        int usedExp = 0;
        
        String[] Res = new String[ResultsSize];
        Res[3] = "-"; Res[4] = "-"; Res[5] = "-"; Res[6] = "-";

        for (InsilicoModelWrapper Wrapper : CurModels) {
        
            InsilicoModelOutput curResult = Wrapper.getResult().get(CurMoleculeIndex);
            
            // If model has not been calculated or error - skip
            if (curResult.getStatus() != InsilicoModelOutput.OUTPUT_OK)
                continue;
            
            // CAESAR muta model
            if (Wrapper.getModel().getClass() == ismMutagenicityCaesar.class) {
                
                Res[3] = curResult.getAssessment();
                
                if (curResult.HasExperimental()) {
                    if (curResult.getExperimental() == 0) {
                        weightNonMutaExp += 1;
                        usedExp++;
                    } 
                    if (curResult.getExperimental() == 1) {
                        weightMutaExp += 1;
                        usedExp++;
                    }
                } else {
                    
                    // non muta
                    if (curResult.getMainResultValue() == 0) {
                        weightNonMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                    
                    // muta or suspect muta
                    if (curResult.getMainResultValue() > 0) {
                        weightMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                }
            }
            
            // ISS muta model
            if (Wrapper.getModel().getClass() == ismMutagenicityBB.class) {
                
                Res[4] = curResult.getAssessment();
                
                if (curResult.HasExperimental()) {
                    if (curResult.getExperimental() == 0) {
                        weightNonMutaExp += 1;
                        usedExp++;
                    } 
                    if (curResult.getExperimental() == 1) {
                        weightMutaExp += 1;
                        usedExp++;
                    }
                } else {
                    
                    // non muta
                    if (curResult.getMainResultValue() == 0) {
                        weightNonMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                    
                    // muta 
                    if (curResult.getMainResultValue() > 0) {
                        weightMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                }
            }
            
            // SARPY muta model
            if (Wrapper.getModel().getClass() == ismMutagenicitySarpy.class) {
                
                Res[5] = curResult.getAssessment();
                
                if (curResult.HasExperimental()) {
                    if (curResult.getExperimental() == 0) {
                        weightNonMutaExp += 1;
                        usedExp++;
                    } 
                    if (curResult.getExperimental() == 1) {
                        weightMutaExp += 1;
                        usedExp++;
                    }
                } else {
                    
                    // non muta and suspect non-muta
                    if (curResult.getMainResultValue() < 1) {
                        weightNonMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                    
                    // muta or suspect muta
                    if (curResult.getMainResultValue() > 0) {
                        weightMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                }
            }
            
            // KNN muta model
            if (Wrapper.getModel().getClass() == ismMutagenicityKnn.class) {
                
                Res[6] = curResult.getAssessment();
                
                if (curResult.HasExperimental()) {
                    if (curResult.getExperimental() == 0) {
                        weightNonMutaExp += 1;
                        usedExp++;
                    } 
                    if (curResult.getExperimental() == 1) {
                        weightMutaExp += 1;
                        usedExp++;
                    }
                } else {
                    
                    // non muta 
                    if (curResult.getMainResultValue() == 0) {
                        weightNonMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                    
                    // muta or suspect muta
                    if (curResult.getMainResultValue() > 0) {
                        weightMuta += GetPredictionWeight(curResult);
                        usedModels++;
                    } 
                }
            }            
            
        }
        
        if ( (usedModels + usedExp) == 0)
            return MODEL_ERROR;
        
        // if at least one exp is retrieved from a model
        // only exp values are used for the consensus
        if (usedExp > 0) {
            usedModels = usedExp;
            weightMuta = weightMutaExp;
            weightNonMuta = weightNonMutaExp;
            CurOutput.setExperimentalBased(true);
        }
        
        // normalize values
        weightMuta = weightMuta / usedModels;
        weightNonMuta = weightNonMuta / usedModels;
        String wM = Format_3D.format(weightMuta);
        String wNM = Format_3D.format(weightNonMuta);
        weightMuta = Double.valueOf(wM);
        weightNonMuta = Double.valueOf(wNM);
        
        CurOutput.setUsedModels(usedModels);
        if (Math.abs(weightMuta - weightNonMuta) < 0.001)  {
            CurOutput.setMainResultValue(1);
            Res[0] = "Mutagenic";
        } else { 
            if (weightMuta > weightNonMuta) {
                CurOutput.setMainResultValue(1);
                Res[0] = "Mutagenic";
            } else {
                CurOutput.setMainResultValue(0);
                Res[0] = "NON-Mutagenic";
            }
        }
        Res[1] = wM;
        Res[2] = wNM;
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;    
    
    }

    @Override
    protected void CalculateAssessment() {
        
        CurOutput.setAssessment(CurOutput.getResults()[0] + " (Consensus score: " + (CurOutput.getMainResultValue()==1?CurOutput.getResults()[1]:CurOutput.getResults()[2]) + ")");
        CurOutput.setAssessmentVerbose("Prediction is " + CurOutput.getResults()[0] + 
                " with a consensus score of " + (CurOutput.getMainResultValue()==1?CurOutput.getResults()[1]:CurOutput.getResults()[2]) + ", based on " +
                CurOutput.getUsedModels() + (CurOutput.isExperimentalBased()? " experimental values." : " models.") );
        
        if (CurOutput.getMainResultValue() == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        if (CurOutput.getMainResultValue() == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        
    }

    @Override
    public ArrayList<InsilicoModel> GetRequiredModels() throws InitFailureException {
        ArrayList<InsilicoModel> models = new ArrayList<>();
        models.add(new ismMutagenicityCaesar());
        models.add(new ismMutagenicityBB());
        models.add(new ismMutagenicitySarpy());
        models.add(new ismMutagenicityKnn());
        return models;
    }
    
}
