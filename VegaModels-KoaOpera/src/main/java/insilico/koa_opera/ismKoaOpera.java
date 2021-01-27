package insilico.koa_opera;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import libpadeldescriptor.MLFERDescriptor;
import libpadeldescriptor.PaDELHBondDonorCountDescriptor;
import opera_adapter.OperaDistance;
import opera_adapter.OperaModel;
import org.openscience.cdk1.qsar.IMolecularDescriptor;
import padeladapter.PadelInterface;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro
 */
public class ismKoaOpera extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_koa_opera.xml";

    private final PadelInterface Padel;
    
    
    private final double[] Scaling_a = {
        0.11481481481481481,
        6.16557037037037  
    };
    
    private final double[] Scaling_s = {
        0.3418772836898127,
        3.476624490000608 
    };
    
    
    public ismKoaOpera() 
            throws InitFailureException {
        super(ModelData);
        
        // Calculator for PaDel descriptors
        Padel = new PadelInterface();
        
        // Define no. of descriptors
        this.DescriptorsSize = 2;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nHBDon (autoscaled)";
        this.DescriptorsNames[1] = "MLFER_L (autoscaled)";
                     
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted KOA [log units]";
        
        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        
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

            //nHBDon
            MD = (IMolecularDescriptor) new PaDELHBondDonorCountDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("nHBDon");
            double[] res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[0] = res[0];
            
            //nHBAcc
            MD = (IMolecularDescriptor) new MLFERDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("MLFER_L");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[1] = res[0];

            
            // Autoscaling of descriptors
            for (int i=0; i<DescriptorsSize; i++)
                Descriptors[i] = (Descriptors[i] - Scaling_a[i]) / Scaling_s[i];
            
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {
        
        // Build model object 
        int Opera_K = 5;
        short Opera_Distance = OperaDistance.DIST_EUCLIDEAN;
        OperaModel OperaKNN = new OperaModel(TS, 2, Opera_K, Opera_Distance);
        
        double OperaResult = 0;
        
        try {
            OperaResult = OperaKNN.Calculate(Descriptors);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }
        
        CurOutput.setMainResultValue(OperaResult);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_4D.format(OperaResult)); // log units
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        // Ad is performed on the K molecules used for the KNN model
        adq.setMoleculesForIndexSize(3);

        // (only retrieve similar molecules if n.a. prediction)
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == Descriptor.MISSING_VALUE) {
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
        // (values are lower than other quantitative models because the data
        //  are expressed in mol NOT mg)
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.75, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.5);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.5);
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

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.75, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        
        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set (mg/L) if available

        String ADItemWarnings =
                ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());
        
        String Result = CurOutput.getResults()[0] + " log units";

        if (CurOutput.getMainResultValue() == Descriptor.MISSING_VALUE) {
            CurOutput.setAssessment("N/A");
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_NA, "N/A"));
        } else {
            switch (CurOutput.getADI().GetAssessmentClass()) {
                case ADIndex.INDEX_LOW:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_LOW, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_LOW, Result, ADItemWarnings));
                    break;
                case ADIndex.INDEX_MEDIUM:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_MEDIUM, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_MEDIUM, Result, ADItemWarnings));
                    break;
                case ADIndex.INDEX_HIGH:
                    CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_HIGH, Result));
                    CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_HIGH, Result));
                    if (!ADItemWarnings.isEmpty())
                        CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() +
                                String.format(MessagesAD.ASSESS_LONG_ADD_ISSUES, ADItemWarnings));
                    break;
            }
        }
        
        // Override assessment if experimental value is available
        if (CurOutput.HasExperimental()) {
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getExperimentalFormatted() + " log units", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getExperimentalFormatted() + " log units"));
        }
        

        // Sets assessment status - no thresholds for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
//        TSK.SetCalculateDescriptors(true);
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this, true, true);
        TSK.SerializeToFile(DatName);        
    }

}
