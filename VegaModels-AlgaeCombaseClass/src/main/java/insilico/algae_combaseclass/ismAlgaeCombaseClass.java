package insilico.algae_combaseclass;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import protoqsar.algae.AlgaeClassificationANN;
import protoqsar.algae.AlgaeDescriptors;
import protoqsar.filters.BiocideFilter;

import java.util.ArrayList;


public class ismAlgaeCombaseClass extends InsilicoModel {

    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(ismAlgaeCombaseClass.class);

    private static final String ModelData = "/data/model_algae_combaseClass.xml";


    public ismAlgaeCombaseClass()
            throws InitFailureException {
        super(ModelData);


        // Define no. of descriptors
        this.DescriptorsSize = 14;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ATS5m";
        this.DescriptorsNames[1] = "B01[C-Cl]";
        this.DescriptorsNames[2] = "B01[C-O]";
        this.DescriptorsNames[3] = "B02[Cl-Cl]";
        this.DescriptorsNames[4] = "B02[N-O]";
        this.DescriptorsNames[5] = "B03[N-O]";
        this.DescriptorsNames[6] = "B09[C-C]";
        this.DescriptorsNames[7] = "B09[O-O]";
        this.DescriptorsNames[8] = "B10[C-O]";
        this.DescriptorsNames[9] = "F01[C-O]";
        this.DescriptorsNames[10] = "F02[C-O]";
        this.DescriptorsNames[11] = "F04[O-O]";
        this.DescriptorsNames[12] = "F10[C-N]";
        this.DescriptorsNames[13] = "X3A";

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted toxicity activity";
        this.ResultsName[1] = "Biocide filter";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
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

            AlgaeDescriptors AlgaeDescEngine = new AlgaeDescriptors();
            AlgaeDescEngine.Calculate(this.CurMolecule);
            double[] algaeDesc = AlgaeDescEngine.GetDescriptorsForAlgaeClassificationModel();

            Descriptors = new double[DescriptorsSize];
            if (DescriptorsSize >= 0) System.arraycopy(algaeDesc, 0, Descriptors, 0, DescriptorsSize);

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }


    @Override
    protected short CalculateModel() {

        AlgaeClassificationANN AlgaeANN = new AlgaeClassificationANN();
        AlgaeANN.RunPrediction(Descriptors);
        int MainResult = AlgaeANN.ResultClass;

        CurOutput.setMainResultValue(MainResult);

        String biocide = BiocideFilter.ApplyFilter(CurMolecule);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // toxicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for algae toxicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = biocide;

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(3);
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
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }


}
