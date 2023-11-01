package insilico.tissueblood_ineris;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.pmml.ModelGenericFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import libpadeldescriptor.AutocorrelationDescriptor;
import libpadeldescriptor.EStateAtomTypeDescriptor;
import libpadeldescriptor.PaDELWeightedPathDescriptor;
import org.dmg.pmml.FieldName;
import org.openscience.cdk1.qsar.IMolecularDescriptor;
import org.openscience.cdk1.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk1.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk1.qsar.descriptors.molecular.XLogPDescriptor;
import padeladapter.PadelInterface;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 *
 * @author User
 */
public class ismTissueBloodIneris extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/tissuebloodineris/model_tissueblood_ineris.xml";
    
    private final double[] MLR_COEFF = {
        0.3393, // AATS8m
        0.26676, //  SssCH2
        -0.43914, // minHBd
        -0.55004, // LipoaffinityIndex
        0.8852, //  XLogP
        0.62931 // intercept       
    };
    
    private final double[] MLR_SCALING_MEAN = {
        92.6744284043791,
        2.49257611594952,
        0.215015334672049,
        5.80376559286982,
        3.18869841269841        
    };
    
    private final double[] MLR_SCALING_STDEV = {
        179.112717847862,
        3.50293070335418,
        0.245964965150894,
        2.44736434663244,
        2.08997105733101        
    };
    
    private final PadelInterface Padel;
    private ModelGenericFromPMML ModelRF;
    
    
    public ismTissueBloodIneris() 
            throws InitFailureException {
        super(ModelData);
        
        // Calculator for PaDel descriptors
        Padel = new PadelInterface();
        
        // Init PMML models
        try {  
            URL src = getClass().getResource("/data/tissuebloodineris/tissueblood_randomforest_model.xml");
            ModelRF = new ModelGenericFromPMML(src.openStream());
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 11;
        this.DescriptorsNames = new String[DescriptorsSize];

        // RF descriptors
        this.DescriptorsNames[0] = "ALogp2";
        this.DescriptorsNames[1] = "ATSC1s";
        this.DescriptorsNames[2] = "BCUTp.1l";
        this.DescriptorsNames[3] = "minHBd";
        this.DescriptorsNames[4] = "XLogP";
        this.DescriptorsNames[5] = "WTPT.5";
        
        // MLR descriptors
        this.DescriptorsNames[6] = "AATS8m";
        this.DescriptorsNames[7] = "SssCH2";
        this.DescriptorsNames[8] = "minHBd";
        this.DescriptorsNames[9] = "LipoaffinityIndex";
        this.DescriptorsNames[10] = "XLogP";
        

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LogKab [log units]";
        
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

            // ALogp2
            MD = new ALOGPDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("ALogp2");
            double[] res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[0] = res[0]; 
            
            // ATSC1s, AATS8m
            MD = new AutocorrelationDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("ATSC1s");
            DescNames.add("AATS8m");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[1] = res[0];
            Descriptors[6] = res[1];

            // BCUTp.1l
            MD = new BCUTDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("BCUTp-1l");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[2] = res[0];
            
            // minHBd, SssCH2, LipoaffinityIndex
            MD = new EStateAtomTypeDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("minHBd");
            DescNames.add("SssCH2");
            DescNames.add("LipoaffinityIndex");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[3] = res[0];
            Descriptors[8] = res[0];
            Descriptors[7] = res[1];
            Descriptors[9] = res[2];            

            // XLogP
            MD = new XLogPDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("XLogP");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[4] = res[0];
            Descriptors[10] = res[0];
            
            // WTPT.5
            MD = new PaDELWeightedPathDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("WTPT-5");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[5] = res[0];
                        
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        //// RF
        
        // no scaling needed on RF descriptors
        
        // Set descriptors and run RF (from idx=0 to idx=5)
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        for (int i=0; i<6; i++)
            argumentsObject.put(this.DescriptorsNames[i], Descriptors[i]);
        
        // Run pmml model
        double PredictionRF;
        try {
            Map<FieldName, ?> OutPred = ModelRF.Evaluate(argumentsObject);
            PredictionRF = (Double)OutPred.get(FieldName.create(("Predicted_LogKab")));
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        
        //// MLR
        
        // Scaling of descriptors (autocentering on training stats)
        double[] mlr_scaled_desc = new double[5];
        for (int i=0; i<5; i++)
            mlr_scaled_desc[i] = (Descriptors[6+i] - MLR_SCALING_MEAN[i]) / MLR_SCALING_STDEV[i];
        
        double PredictionMLR = 0;
        for (int i=0; i<5; i++)
            PredictionMLR += mlr_scaled_desc[i] * MLR_COEFF[i];
        PredictionMLR +=  MLR_COEFF[5];
        
        
        //// Final prediction - from RF model
        
        double Prediction = PredictionRF;
        
        
        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_4D.format(Prediction));

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
