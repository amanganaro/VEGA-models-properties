package insilico.bcf_knn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.knn.insilicoKnnQuantitative;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;


/**
 *
 * @author User
 */
public class ismBCFKnn extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_bcf_knn.xml";

    private final insilicoKnnQuantitative KNN;
    private insilicoKnnPrediction KnnPrediction;
    
    
    public ismBCFKnn() 
            throws InitFailureException {
        super(ModelData);
        
        // Build model object
        KNN =  new insilicoKnnQuantitative();
        KNN.setNeighboursNumber(4);
        KNN.setMinSimilarity(0.7);
        KNN.setMinSimilarityForSingleResult(0.75);
        KNN.setEnhanceWeightFactor(3);
        KNN.setExperimentalRange(3.5);
        KNN.setUseExperimentalRange(true);
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted BCF [log(L/kg)]";
        this.ResultsName[1] = "Molecules used for prediction";
        
        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {
        
        try {
            KnnPrediction = KNN.Calculate(CurMolecule, TS);
        } catch (GenericFailureException ex) {
            return MODEL_ERROR;
        }
        
        CurOutput.setMainResultValue(KnnPrediction.getPrediction());

        String[] Res = new String[ResultsSize];
        if (KnnPrediction.getPrediction() == Descriptor.MISSING_VALUE)
            Res[0] = "-";
        else
            Res[0] = String.valueOf(Format_2D.format(KnnPrediction.getPrediction())); // BCF
        Res[1] = String.valueOf(KnnPrediction.getNeighbours().size());
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        // Ad is performed on the K molecules used for the KNN model
        adq.setMoleculesForIndexSize(KnnPrediction.getNeighbours().size());

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
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.75, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.2, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.2, 0.6);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.2, 0.6);
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
        
        double BCF_threshold_green = 2.7;
        double BCF_threshold_red = 3.3;
    
        // Sets assessment message
        if (CurOutput.getMainResultValue() == Descriptor.MISSING_VALUE) {
            
            CurOutput.setAssessment("N/A");
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_NA, "N/A"));
            
        } else {
    
            ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0], "log(L/kg)");

        }        
        
        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == Descriptor.MISSING_VALUE)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else {
            if (Val < BCF_threshold_green)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
            else if (Val < BCF_threshold_red)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
            else
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        }
    }
    
    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this, true, false);
        TSK.SerializeToFile(DatName);        
    }
    
}
