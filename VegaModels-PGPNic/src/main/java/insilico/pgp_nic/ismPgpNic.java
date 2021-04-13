package insilico.pgp_nic;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.ADIndexACF;
import insilico.core.ad.item.ADIndexADI;
import insilico.core.ad.item.ADIndexAccuracy;
import insilico.core.ad.item.ADIndexConcordance;
import insilico.core.ad.item.ADIndexRange;
import insilico.core.ad.item.ADIndexSimilarity;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.pgp_nic.descriptors.modelCPANN;
import insilico.pgp_nic.descriptors.modelPgpDescriptors;
import lombok.extern.slf4j.Slf4j;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

/**
 *
 * @author User
 */
@Slf4j
public class ismPgpNic extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_pgp_nic.xml";
    
    private final modelCPANN CPANN;
    private final DecimalFormat df2;
    private final DecimalFormat df3;
    
    public ismPgpNic() 
            throws InitFailureException {
        super(ModelData);
        
        DecimalFormatSymbols InternationalSymbols =
            new DecimalFormatSymbols();
        InternationalSymbols.setDecimalSeparator('.');         
        df2 = new DecimalFormat("0.##", InternationalSymbols);
        df3 = new DecimalFormat("0.###", InternationalSymbols);
        
        // Builds model objects
        CPANN = new modelCPANN();
        
        // Define no. of descriptors
        this.DescriptorsSize = 26;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "H%";
        this.DescriptorsNames[1] = "nR07";
        this.DescriptorsNames[2] = "D/Dtr11";
        this.DescriptorsNames[3] = "MWC01";
        this.DescriptorsNames[4] = "X2A";
        this.DescriptorsNames[5] = "SIC3";
        this.DescriptorsNames[6] = "VE1sign_B(s)";
        this.DescriptorsNames[7] = "ATSC7m";
        this.DescriptorsNames[8] = "MATS6v";
        this.DescriptorsNames[9] = "GATS4s";
        this.DescriptorsNames[10] = "P_VSA_LogP_3";
        this.DescriptorsNames[11] = "P_VSA_ppp_D";
        this.DescriptorsNames[12] = "nRCOOR";
        this.DescriptorsNames[13] = "nArCONHR";
        this.DescriptorsNames[14] = "nArCO";
        this.DescriptorsNames[15] = "H-048";
        this.DescriptorsNames[16] = "SdsCH";
        this.DescriptorsNames[17] = "CATS2D_01_DN";
        this.DescriptorsNames[18] = "CATS2D_05_PP";
        this.DescriptorsNames[19] = "CATS2D_02_PL";
        this.DescriptorsNames[20] = "B07[O-F]";
        this.DescriptorsNames[21] = "F01[C-C]";
        this.DescriptorsNames[22] = "F02[C-O]";
        this.DescriptorsNames[23] = "F04[C-P]";
        this.DescriptorsNames[24] = "F04[C-Br]";
        this.DescriptorsNames[25] = "F07[O-F]";        
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted P-gp class";
        this.ResultsName[1] = "Probability for Inhibitor class";
        this.ResultsName[2] = "Probability for Substrate class";
        this.ResultsName[3] = "Probability for Non Active class";
        this.ResultsName[4] = "Euclidean Distance from the central neuron";
        
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
            
            // custom descriptors
            modelPgpDescriptors pgp = new modelPgpDescriptors();
            pgp.Calculate(CurMolecule, null, null);
            
            Descriptors = new double[DescriptorsSize];

            for (int i=0; i<pgp.Descriptors.length; i++)
                this.Descriptors[i] = pgp.Descriptors[i];
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        double[] PredClass = CPANN.CalculatePrediction(this.Descriptors);
        
        // three classes available, final prediction assigned to the class
        // with prob > 0.5 or not predicted if no class has prob > 0.5
        int Prediction = -1;
        for (int i=0; i<3; i++)
            if (PredClass[i] > 0.5) {
                Prediction = i;
                break;
            }

        CurOutput.setMainResultValue(Prediction);
        
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction); // P-gp classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for P-gp value " + Prediction);
            Res[0] = Integer.toString(Prediction);
        }
        Res[1] = df2.format(PredClass[0]);
        Res[2] = df2.format(PredClass[1]);
        Res[3] = df2.format(PredClass[2]);
        Res[4] = df3.format(PredClass[5]);
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS, false);
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
        if ( (Val == 0) || (Val == 1))
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val == 2)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        

    }
}
