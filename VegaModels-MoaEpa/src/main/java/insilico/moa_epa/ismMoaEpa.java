package insilico.moa_epa;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.moa_epa.descriptors.EmbeddedDescriptors;
import insilico.moa_epa.descriptors.MOAToxAdditionalDescriptors;
import insilico.moa_epa.descriptors.MOAToxMultipleModels;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class ismMoaEpa extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_moa_epa.xml";
    
    private final static double SCORE_MIN_THRESHOLD = 0.5;
    private final MOAToxAdditionalDescriptors AdditionalDescriptors;
    
    public ismMoaEpa() 
            throws InitFailureException {
        super(ModelData);
        
        // Init additional descriptors
        AdditionalDescriptors = new MOAToxAdditionalDescriptors();

        // Define no. of descriptors
        this.DescriptorsSize = MOAToxMultipleModels.DescriptorsSize;
        this.DescriptorsNames = new String[DescriptorsSize];
        DescriptorsNames[0] = "SdsssP_acnt";
        DescriptorsNames[1] = "ic"; 
        DescriptorsNames[2] = "MDEC22";
        DescriptorsNames[3] = "MATS1m";
        DescriptorsNames[4] = "MATS2e";
        DescriptorsNames[5] = "GATS2e";
        DescriptorsNames[6] = "-O- [2 phosphorus attach]";
        DescriptorsNames[7] = "-OH [phosphorus attach]";
        DescriptorsNames[8] = "-NH- [aliphatic attach]";
        DescriptorsNames[9] = "-NH- [aromatic attach]";
        DescriptorsNames[10] = "-C(=O)O- [nitrogen attach]";
        DescriptorsNames[11] = "xch10"; // manca X chain 10
        DescriptorsNames[12] = "MDEC13";
        DescriptorsNames[13] = "MDEC23";
        DescriptorsNames[14] = "BELe1";
        DescriptorsNames[15] = "nR10";
        DescriptorsNames[16] = "GATS6v";
        DescriptorsNames[17] = "-CH< [aromatic attach]";
        DescriptorsNames[18] = ">C= [aromatic attach]";
        DescriptorsNames[19] = "C=O(ketone, aliphatic attach)";
        DescriptorsNames[20] = "-C(=O)- [aromatic attach]";
        DescriptorsNames[21] = "-C(=O)O- [cyclic]";
        DescriptorsNames[22] = "Gmin";
        DescriptorsNames[23] = "-C#N [aliphatic nitrogen attach]";
        DescriptorsNames[24] = "-NO2 [nitrogen attach]";
        DescriptorsNames[25] = "SdssNp";
        DescriptorsNames[26] = "SdCH2_acnt";
        DescriptorsNames[27] = "StsC_acnt";
        DescriptorsNames[28] = "SdsssP_acnt";
        DescriptorsNames[29] = "MDEC33";
        DescriptorsNames[30] = "AMW";
        DescriptorsNames[31] = "SRW07";
        DescriptorsNames[32] = "-CH< [aromatic attach]";
        DescriptorsNames[33] = "-NH- [nitrogen attach]";
        DescriptorsNames[34] = "-CHO [aliphatic attach]";
        DescriptorsNames[35] = "-C(=O)O- [nitrogen attach]";   
        DescriptorsNames[36] = "SdssC";
        DescriptorsNames[37] = "StN";
        DescriptorsNames[38] = "ic"; 
        DescriptorsNames[39] = "icycem";
        DescriptorsNames[40] = "MDEC34";
        DescriptorsNames[41] = "nX";
        DescriptorsNames[42] = "nR05";
        DescriptorsNames[43] = "nR06";
        DescriptorsNames[44] = "nR10";
        DescriptorsNames[45] = "SRW07";
        DescriptorsNames[46] = "-C(=O)O- [aliphatic attach]";            
        DescriptorsNames[47] = "StsC";
        DescriptorsNames[48] = "SsOm";
        DescriptorsNames[49] = "SdCH2_acnt";
        DescriptorsNames[50] = "StCH_acnt";
        DescriptorsNames[51] = "BEHp4";
        DescriptorsNames[52] = "-NH- [nitrogen attach]";
        DescriptorsNames[53] = "-S- [sulfur attach]";
        DescriptorsNames[54] = "-C(=O)- [aromatic attach]";
        DescriptorsNames[55] = "-CHO [aliphatic attach]";
        DescriptorsNames[56] = "-CHO [aromatic attach]";
        DescriptorsNames[57] = "CH2=CHC(=O)O-";
        DescriptorsNames[58] = "SdssNp";
        DescriptorsNames[59] = "Hmin";
        DescriptorsNames[60] = "Qsv"; 
        DescriptorsNames[61] = "ivdem"; 
        DescriptorsNames[62] = "Ms"; 
        DescriptorsNames[63] = "ATS4m";
        DescriptorsNames[64] = "GATS8m";
        DescriptorsNames[65] = "-Br [aromatic attach]";
        DescriptorsNames[66] = "-Cl [aromatic attach]";
        DescriptorsNames[67] = "-I [aromatic attach]";
        DescriptorsNames[68] = "-OH [aromatic attach]";
        
        // Defines results
        this.ResultsSize = 1 + MOAToxMultipleModels.MOAModelsNumber;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted MOA toxicity";
        for (int i=0; i<MOAToxMultipleModels.MOAModelsNumber; i++)
            this.ResultsName[1 + i] = "Score for MOA: " + MOAToxMultipleModels.MOANames[i];
        
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
            
            // calculate additional (custom) descriptors
            AdditionalDescriptors.Calculate(CurMolecule);
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];
            
            Descriptors[0] = AdditionalDescriptors.GetByName("SdsssP_acnt");
            Descriptors[1] = AdditionalDescriptors.GetByName("ic");
            Descriptors[2] = embeddedDescriptors.getMDE_C_22();
            Descriptors[3] = embeddedDescriptors.getMATS1m();
            Descriptors[4] = embeddedDescriptors.getMATS2e();
            Descriptors[5] = embeddedDescriptors.getGATS2e();
            Descriptors[6] = AdditionalDescriptors.GetByName("-O- [2 phosphorus attach]");
            Descriptors[7] = AdditionalDescriptors.GetByName("-OH [phosphorus attach]");
            Descriptors[8] = AdditionalDescriptors.GetByName("-NH- [aliphatic attach]");
            Descriptors[9] = AdditionalDescriptors.GetByName("-NH- [aromatic attach]");
            Descriptors[10] = AdditionalDescriptors.GetByName("-C(=O)O- [nitrogen attach]");
       
            Descriptors[11] = AdditionalDescriptors.GetByName("Xch10");
            Descriptors[12] = embeddedDescriptors.getMDE_C_13();
            Descriptors[13] = embeddedDescriptors.getMDE_C_23();
            Descriptors[14] = (-1.0) * embeddedDescriptors.getBEL1e();
            Descriptors[15] = embeddedDescriptors.getNR10();
            Descriptors[16] = embeddedDescriptors.getGATS6v();
            Descriptors[17] = AdditionalDescriptors.GetByName("-CH< [aromatic attach]");
            Descriptors[18] = AdditionalDescriptors.GetByName(">C= [aromatic attach]");
            Descriptors[19] = AdditionalDescriptors.GetByName("C=O(ketone, aliphatic attach)");
            Descriptors[20] = AdditionalDescriptors.GetByName("-C(=O)- [aromatic attach]");
            Descriptors[21] = AdditionalDescriptors.GetByName("-C(=O)O- [cyclic]");
        
            Descriptors[22] = embeddedDescriptors.getGmin();
            Descriptors[23] = AdditionalDescriptors.GetByName("-C#N [aliphatic nitrogen attach]");
            Descriptors[24] = AdditionalDescriptors.GetByName("-NO2 [nitrogen attach]");
            
            Descriptors[25] = embeddedDescriptors.getSdssNp();
            Descriptors[26] = AdditionalDescriptors.GetByName("SdCH2_acnt");
            Descriptors[27] = AdditionalDescriptors.GetByName("StsC_acnt");
            Descriptors[28] = AdditionalDescriptors.GetByName("SdsssP_acnt");
            Descriptors[29] = embeddedDescriptors.getMDE_C_33();
            Descriptors[30] = AdditionalDescriptors.GetByName("AMW");
            Descriptors[31] = embeddedDescriptors.getSRW7();
            Descriptors[32] = AdditionalDescriptors.GetByName("-CH< [aromatic attach]");
            Descriptors[33] = AdditionalDescriptors.GetByName("-NH- [nitrogen attach]");
            Descriptors[34] = AdditionalDescriptors.GetByName("-CHO [aliphatic attach]");
            Descriptors[35] = AdditionalDescriptors.GetByName("-C(=O)O- [nitrogen attach]");
            
            Descriptors[36] = embeddedDescriptors.getSdssC();
            Descriptors[37] = embeddedDescriptors.getStN();
            Descriptors[38] = AdditionalDescriptors.GetByName("ic");
            Descriptors[39] = AdditionalDescriptors.GetByName("icycem");
            Descriptors[40] = embeddedDescriptors.getMDE_C_34();
            Descriptors[41] = embeddedDescriptors.getNX();
            Descriptors[42] = embeddedDescriptors.getNR5();
            Descriptors[43] = embeddedDescriptors.getNR6();
            Descriptors[44] = embeddedDescriptors.getNR10();
            Descriptors[45] = embeddedDescriptors.getSRW7();
            Descriptors[46] = AdditionalDescriptors.GetByName("-C(=O)O- [aliphatic attach]");

            Descriptors[47] = embeddedDescriptors.getStsC();
            Descriptors[48] = embeddedDescriptors.getSsOm();
            Descriptors[49] = AdditionalDescriptors.GetByName("SdCH2_acnt");
            Descriptors[50] = AdditionalDescriptors.GetByName("StCH_acnt");
            Descriptors[51] = embeddedDescriptors.getBEH4p();
            Descriptors[52] = AdditionalDescriptors.GetByName("-NH- [nitrogen attach]");
            Descriptors[53] = AdditionalDescriptors.GetByName("-S- [sulfur attach]");
            Descriptors[54] = AdditionalDescriptors.GetByName("-C(=O)- [aromatic attach]");
            Descriptors[55] = AdditionalDescriptors.GetByName("-CHO [aliphatic attach]");
            Descriptors[56] = AdditionalDescriptors.GetByName("-CHO [aromatic attach]");
            Descriptors[57] = AdditionalDescriptors.GetByName("CH2=CHC(=O)O-");
            
            Descriptors[58] = embeddedDescriptors.getSdssNp();
            Descriptors[59] = embeddedDescriptors.getHmin();
            // in TEST, when the descriptor is not calculated it is set to 0
            Descriptors[59] = (Descriptors[59] == -999 ? 0.0 : Descriptors[59]);
            Descriptors[60] = embeddedDescriptors.getQsv();
            Descriptors[61] = embeddedDescriptors.getIVDEM();
            Descriptors[62] = embeddedDescriptors.getMs();
            Descriptors[63] = embeddedDescriptors.getATS4m();
            Descriptors[64] = embeddedDescriptors.getGATS8m();
            Descriptors[65] = AdditionalDescriptors.GetByName("-Br [aromatic attach]");
            Descriptors[66] = AdditionalDescriptors.GetByName("-Cl [aromatic attach]");
            Descriptors[67] = AdditionalDescriptors.GetByName("-I [aromatic attach]");
            Descriptors[68] = AdditionalDescriptors.GetByName("-OH [aromatic attach]");
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        int MainResult = -1;
        
        // Calculate scores for all MOA classes
        double[] Scores = MOAToxMultipleModels.CalculateScores(Descriptors);
        
        // Find best matching class
        int MaxScoreIndex = 0;
        for (int i=0; i<Scores.length; i++) 
            if (Scores[MaxScoreIndex] < Scores[i])
                MaxScoreIndex = i;
        
        // Decide final MOA class
        if (Scores[MaxScoreIndex] <= SCORE_MIN_THRESHOLD)
            MainResult = -1;
        else
            MainResult = MaxScoreIndex;
        
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // MOA class
        } catch (Throwable ex) {
            log.warn("Unable to find label for carcinogenicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        for (int i=0; i<MOAToxMultipleModels.MOAModelsNumber; i++)
            Res[1 + i] = Format_4D.format(Scores[i]);
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS, false);
        // Model with multiple classes, no value mapping is used
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

        // No assessment color for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        
    }
    
}
