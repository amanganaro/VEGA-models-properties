package insilico.totalhl_qsarins;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import libpadeldescriptor.AutocorrelationDescriptor;
import libpadeldescriptor.EStateAtomTypeDescriptor;
import libpadeldescriptor.HalogenCountDescriptor;
import libpadeldescriptor.TopologicalChargeDescriptor;
import org.openscience.cdk1.qsar.IMolecularDescriptor;
import org.openscience.cdk1.qsar.descriptors.molecular.TPSADescriptor;
import padeladapter.PadelInterface;

import java.util.ArrayList;


/**
 *
 * @author User
 */
public class ismTotalHLQsarins extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_totalhl_qsarins.xml";
    
    private final double[] MLR_COEFF_FULLMODEL = {
        0.2105, //nF
        - 0.2927, //minHsOH
        - 0.0495, //gmax
        0.3510, //SsCl
        0.5905, //AATS7p
        - 0.0042, //TopoPSA
        0.0686, //GGI1
        - 0.4298, //minsCl
        0.6577        
    };
    
    private final double[] MLR_COEFF_SPLIT ={
        + 0.2385, // nF 
        - 0.27, // minsHsOH
        - 0.0484, // gmax 
        + 0.3299, // ScCl 
        + 0.6018, // AATS7p 
        - 0.0043, // TopoPSA 
        + 0.0778, // GGI1 
        - 0.2404, // minsCl 
        0.5683        
    };
    
    private final PadelInterface Padel;
    
    
    public ismTotalHLQsarins() 
            throws InitFailureException {
        super(ModelData);
        
        // Calculator for PaDel descriptors
        Padel = new PadelInterface();
        
        
        // Define no. of descriptors
        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nF";
        this.DescriptorsNames[1] = "minHsOH";
        this.DescriptorsNames[2] = "gmax";
        this.DescriptorsNames[3] = "SsCl";
        this.DescriptorsNames[4] = "AATS7p";
        this.DescriptorsNames[5] = "TopoPSA";
        this.DescriptorsNames[6] = "GGI1";
        this.DescriptorsNames[7] = "minsCl";
        

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LogHLt [log units]";
        this.ResultsName[1] = "Predicted total half-life [hours]";
        
        // Define AD items
        this.ADItemsName = new String[6];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexRange().GetIndexName();
        this.ADItemsName[5] = new ADIndexACF().GetIndexName();

    }
    
    
    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
        
        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        
        return blocks;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
            
            IMolecularDescriptor MD;
            ArrayList<String> DescNames;

            Padel.SetSMILES(CurMolecule.GetSMILES());

            // nF
            MD = new HalogenCountDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("nF");
            double[] res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[0] = res[0]; 
            
            // minHsOH, gmax, SsCl, minsCl
            MD = new EStateAtomTypeDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("minHsOH");
            DescNames.add("gmax");
            DescNames.add("SsCl");
            DescNames.add("minsCl");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[1] = res[0];
            Descriptors[2] = res[1];
            Descriptors[3] = res[2];             
            Descriptors[7] = res[3];             
            
            // AATS7p
            MD = new AutocorrelationDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("AATS7p");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[4] = res[0];
            
            // TopoPSA
            MD = new TPSADescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("TopoPSA");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[5] = res[0];
            
            // GGI1
            MD = new TopologicalChargeDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("GGI1");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[6] = res[0];
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        // Linear regression (Split Model)
        double PredictionSplit = 0;
        for (int i=0; i<8; i++)
            PredictionSplit += Descriptors[i] * MLR_COEFF_SPLIT[i];
        PredictionSplit +=  MLR_COEFF_SPLIT[8];
        
        // Linear regression (Full Model)
        double PredictionFull = 0;
        for (int i=0; i<8; i++)
            PredictionFull += Descriptors[i] * MLR_COEFF_FULLMODEL[i];
        PredictionFull +=  MLR_COEFF_FULLMODEL[8];
        
        // Prediction
        // Final prediction is just the full model
        double Prediction = PredictionFull;
        
        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_3D.format(Prediction));
        double ConvertedValue = Math.pow(10, Prediction);
        if (ConvertedValue > 10)
            Res[1] = String.valueOf(Format_2D.format(ConvertedValue));
        else if (ConvertedValue > 0.01)
            Res[1] = String.valueOf(Format_3D.format(ConvertedValue));
        else
            Res[1] = String.valueOf(Format_4D.format(ConvertedValue));

        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.5, 0.8);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.5, 0.8);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.5, 0.8);
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

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);
        
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status - no thresholds for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
