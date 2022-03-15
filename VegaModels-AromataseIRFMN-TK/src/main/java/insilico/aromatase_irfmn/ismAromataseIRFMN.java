package insilico.aromatase_irfmn;

import insilico.aromatase_irfmn.descriptors.AromataseDescriptors;
import insilico.aromatase_irfmn.descriptors.AromatasePVSA;
import insilico.aromatase_irfmn.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelGenericFromPMML;

import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;
import org.dmg.pmml.FieldName;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author User
 */
@Slf4j
public class ismAromataseIRFMN extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_aromatase_irfmn.xml";
    
    private static final double[] NormalizationStd = {0.711769967075992,0.290475667889515,0.115044003629162,2.32014205188935,2.67041727590661,1.60271810184668,39.4984891454959,2.79017398841024,6.21542806393645,0.824916629012428,1.37409213162412,0.968245511028326,0.413271775530228,0.940353551341537,1.95524050568146,0.584625806488982,3.54631053724478,5.77201335431423};
    private static final double[] NormalizationMean = {3.74823674096849,3.54761913912375,0.236886241352806,4.42416448885473,4.18082013835511,2.10361568024596,56.8000303612606,4.04461606456572,7.87367063797079,2.56390392006149,0.7273904688701,9.01345157571099,0.162567255956956,5.28766372021522,1.05276441199078,0.194081475787856,3.35395849346656,34.6697601844735};

    private final AromatasePVSA pvsa;

    private ModelGenericFromPMML Model;
    
    
    public ismAromataseIRFMN() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML model
        try {  
            URL src = getClass().getResource("/data/rangerRF_aromatase_irfmn.pmml");
            Model = new ModelGenericFromPMML(src.openStream());
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // custom descriptors block
        pvsa = new AromatasePVSA();
        pvsa.setBoolProperty(AromatasePVSA.PARAMETER_WEIGHT_MR, true);
        
        // Define no. of descriptors
        this.DescriptorsSize = 18;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "SpMax2_Bh.m.";
        this.DescriptorsNames[1] = "SpMax2_Bh.p.";
        this.DescriptorsNames[2] = "SpMaxA_AEA.dm.";
        this.DescriptorsNames[3] = "piPC08";
        this.DescriptorsNames[4] = "piPC09";
        this.DescriptorsNames[5] = "MLOGP";
        this.DescriptorsNames[6] = "P_VSA_MR_6";
        this.DescriptorsNames[7] = "X5sol";
        this.DescriptorsNames[8] = "ATSC4p";
        this.DescriptorsNames[9] = "Eig03_EA.bo.";
        this.DescriptorsNames[10] = "SsssN";
        this.DescriptorsNames[11] = "SM6_B.m.";
        this.DescriptorsNames[12] = "nRNR2";
        this.DescriptorsNames[13] = "ATS5s";
        this.DescriptorsNames[14] = "SssNH";
        this.DescriptorsNames[15] = "H.049";
        this.DescriptorsNames[16] = "F02.C.N.";
        this.DescriptorsNames[17] = "SM15_EA.ed.";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Aromatase activity";
        this.ResultsName[1] = "Probability(Active Agonist)";
        this.ResultsName[2] = "Probability(Active Antagonist)";
        this.ResultsName[3] = "Probability(Inactive)";
        
        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

    }
    

    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {

            pvsa.Calculate(CurMolecule);

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.BEH2m;
            Descriptors[1] = embeddedDescriptors.BEH2p;
            Descriptors[2] = AromataseDescriptors.Calculate_spMaxAEAdm(CurMolecule);
            Descriptors[3] = embeddedDescriptors.piPC08;
            Descriptors[4] = embeddedDescriptors.piPC09;
            Descriptors[5] = embeddedDescriptors.MLogP;
            Descriptors[6] = pvsa.GetByName("P_VSA_mr_5").getValue();
            Descriptors[7] = embeddedDescriptors.X5sol;
            Descriptors[8] = embeddedDescriptors.ATSC4p;
            Descriptors[9] = embeddedDescriptors.EEig3bo;
            Descriptors[10] = embeddedDescriptors.SsssN;
            Descriptors[11] = AromataseDescriptors.Calculate_SM6_Bm(CurMolecule);
            Descriptors[12] = embeddedDescriptors.nRNR2;
            Descriptors[13] = embeddedDescriptors.ATS5s;
            Descriptors[14] = embeddedDescriptors.SssNH;
            Descriptors[15] = embeddedDescriptors.H049;
            Descriptors[16] = embeddedDescriptors.F02C_N;
            Descriptors[17] = AromataseDescriptors.Calculate_SM15_EAed(CurMolecule);

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        // Scale input descriptors
        double[] ScaledDescriptors = new double[this.DescriptorsSize];
        for (int i=0; i<DescriptorsSize; i++)
            ScaledDescriptors[i] = (Descriptors[i] - NormalizationMean[i]) / NormalizationStd[i];
        
        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        for (int i=0; i<DescriptorsSize; i++)
            argumentsObject.put(this.DescriptorsNames[i], ScaledDescriptors[i]);
        
        // Run pmml model
        double pInactive, pAgonist, pAntagonist;
        try {
            Map<FieldName, ?> oo = Model.Evaluate(argumentsObject);
            pInactive = (Double)oo.get(FieldName.create("probability(inactive)"));
            pAgonist = (Double)oo.get(FieldName.create("probability(active.agonist)"));
            pAntagonist = (Double)oo.get(FieldName.create("probability(active.antagonist)"));
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        // Decide final class
        double Prediction;
        if ( (pInactive > pAgonist) && (pInactive > pAntagonist) )
            Prediction = 0;
        else if ( (pAgonist > pInactive) && (pAgonist > pAntagonist) )
            Prediction = 1;
        else if ( (pAntagonist > pInactive) && (pAntagonist > pAgonist) )
            Prediction = 2;
        else
            Prediction = -1;

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction); // Carcinogenicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for aromatase activity value " + Prediction);
            Res[0] = String.valueOf(Prediction);
        }
        Res[1] = Format_3D.format(pAgonist); 
        Res[2] = Format_3D.format(pAntagonist); 
        Res[3] = Format_3D.format(pInactive);
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
        adq.AddMappingToNegativeValue(-1);
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
        else if ( (Val == 1) || (Val == 2) )
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else 
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        
    }
    
}
