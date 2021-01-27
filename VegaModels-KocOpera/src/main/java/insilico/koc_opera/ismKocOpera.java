package insilico.koc_opera;

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
import libpadeldescriptor.*;
import opera_adapter.OperaDistance;
import opera_adapter.OperaModel;
import org.openscience.cdk1.qsar.IMolecularDescriptor;
import org.openscience.cdk1.qsar.descriptors.molecular.VABCDescriptor;
import padeladapter.PadelInterface;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro
 */
public class ismKocOpera extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_koc_opera.xml";

    private final PadelInterface Padel;
    
    private final double[] Scaling_a = {
        3.0307891769547357,
        2.270233196159122,
        63.35543537555478,
        5.318413700144006,
        24.810699588477366,
        285.3113854595336,
        196.9813707465152,
        0.31039025557117533,
        7.79835390946502,
        129.53696094024463,
        33.204389574759944,
        7229.489830705381    
    };
    
    private final double[] Scaling_s = {
        1.715441001163294,
        1.9985448210107284,
        49.16433384612706,
        2.251988911908191,
        15.0021919892538,
        336.0057608323125,
        73.77056520895081,
        0.7673630941068451,
        3.1249962691607083,
        63.67950909708482,
        40.91321131158057,
        2960.998453040209  
    };
    
    
    public ismKocOpera() 
            throws InitFailureException {
        super(ModelData);
        
        // Calculator for PaDel descriptors
        Padel = new PadelInterface();
        
        // Define no. of descriptors
        this.DescriptorsSize = 12;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "CrippenLogP (autoscaled)";
        this.DescriptorsNames[1] = "nHBAcc (autoscaled)";
        this.DescriptorsNames[2] = "SpDiam_Dzi (autoscaled)";
        this.DescriptorsNames[3] = "VP-1 (autoscaled)";
        this.DescriptorsNames[4] = "MPC3 (autoscaled)";
        this.DescriptorsNames[5] = "TPC (autoscaled)";
        this.DescriptorsNames[6] = "VABC (autoscaled)";
        this.DescriptorsNames[7] = "SsssN (autoscaled)";
        this.DescriptorsNames[8] = "topoDiameter (autoscaled)";
        this.DescriptorsNames[9] = "AATS1m (autoscaled)";
        this.DescriptorsNames[10] = "MPC6 (autoscaled)";
        this.DescriptorsNames[11] = "ATS1v (autoscaled)";
                     
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted KOC [log(L/Kg)]";
        
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

            //CrippenLogP
            MD = new CrippenDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("CrippenLogP");
            double[] res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[0] = res[0];
            
            //nHBAcc
            MD = new PaDELHBondAcceptorCountDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("nHBAcc");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[1] = res[0];

            //SpDiam_Dzi
            MD = new BaryszMatrixDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("SpDiam_Dzi");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[2] = res[0];

            //VP_1
            MD = new PaDELChiPathDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("VP-1");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[3] = res[0];

            //MPC3
            //TPC
            MD = new PathCountDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("MPC3");
            DescNames.add("TPC");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[4] = res[0];
            Descriptors[5] = res[1];

            //VABC
            MD = new VABCDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("VABC");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[6] = res[0];

            //SsssN
            MD = new EStateAtomTypeDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("SsssN");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[7] = res[0];

            //topoDiameter
            MD = new TopologicalDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("topoDiameter");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[8] = res[0];

            //AATS1m
            MD = new AutocorrelationDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("AATS1m");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[9] = res[0];

            //MPC6
            MD = new PathCountDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("MPC6");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[10] = res[0];

            //ATS1v 
            MD = new AutocorrelationDescriptor();
            DescNames = new ArrayList<>();
            DescNames.add("ATS1v");
            res = Padel.CalculateDescriptors(MD, DescNames);
            Descriptors[11] = res[0];

            
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
        OperaModel OperaKNN = new OperaModel(TS, 12, Opera_K, Opera_Distance);
        
        double OperaResult = 0;
        
        try {
            OperaResult = OperaKNN.Calculate(Descriptors);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }
        
        CurOutput.setMainResultValue(OperaResult);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_4D.format(OperaResult)); // log(L/Kg)
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
        
        String Result = CurOutput.getResults()[0] + " log(L/Kg)";

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
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getExperimentalFormatted() + " log(L/Kg)", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getExperimentalFormatted() + " log(L/Kg)"));
        }
        

        // Sets assessment status - no thresholds for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this, true,true);
        TSK.SerializeToFile(DatName);        
    }

}
