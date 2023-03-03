package insilico.carcinogenicity_caesar;

import insilico.carcinogenicity_caesar.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
//import insilico.core.tools.ModelUtilities;
//import insilico.core.tools.logger.InsilicoLogger;
//


/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class ismCarcinogenicityCaesar extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismCarcinogenicityCaesar.class);


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_carc_caesar.xml";

    private static final double CARC_THRESHOLD = 0.5;
    private final modelCPANN CPANN;
    private final ADNeuronsChecker NeuronChecker;
    private double[] ANNResults;


    public ismCarcinogenicityCaesar()
            throws InitFailureException {
        super(ModelData);

        // Builds model objects
        CPANN = new modelCPANN();
        NeuronChecker = new ADNeuronsChecker(35, 35, "/data/carccaesar_neurons.txt");

        // Define no. of descriptors
        this.DescriptorsSize = 12;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "PW5";
        this.DescriptorsNames[1] = "D/Dr6";
        this.DescriptorsNames[2] = "MATS2p";
        this.DescriptorsNames[3] = "EEig10x";
        this.DescriptorsNames[4] = "ESpm11x";
        this.DescriptorsNames[5] = "ESpm9d";
        this.DescriptorsNames[6] = "GGI2";
        this.DescriptorsNames[7] = "JGI6";
        this.DescriptorsNames[8] = "nRNNOx";
        this.DescriptorsNames[9] = "nPO4";
        this.DescriptorsNames[10] = "N_067";
        this.DescriptorsNames[11] = "N_078";

        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Carcinogen activity";
        this.ResultsName[1] = "P(Carcinogen)";
        this.ResultsName[2] = "P(NON-Carcinogen)";

        // Define AD items
        this.ADItemsName = new String[7];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();
        this.ADItemsName[5] = new ADIndexNNReliability().GetIndexName();
        this.ADItemsName[6] = new ADIndexNNNeuronsConcordance().GetIndexName();

    }





    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();
            embeddedDescriptors.CalculateAllDescriptors(CurMolecule);


            Descriptors[0] = embeddedDescriptors.PW5;
            Descriptors[1] = embeddedDescriptors.DDr6;
            Descriptors[2] = embeddedDescriptors.MATS2p;
            Descriptors[3] = embeddedDescriptors.EEig10ed;
            Descriptors[4] = embeddedDescriptors.ESpm11ed;
            Descriptors[5] = embeddedDescriptors.ESpm9dm;
            Descriptors[6] = embeddedDescriptors.GGI2;
            Descriptors[7] = embeddedDescriptors.JGI6;


            // Functional groups are calculated here
            CarcinogenicityGroups Groups = new CarcinogenicityGroups(CurMolecule);
            Descriptors[8] = Groups.getnRNNOx();
            Descriptors[9] = Groups.getnPO4();
            Descriptors[10] = Groups.getN_067();
            Descriptors[11] = Groups.getN_078();

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }



    @Override
    protected short CalculateModel() {

        ANNResults = CPANN.CalculatePrediction(Descriptors);

        int MainResult;

        if (ANNResults[0] > CARC_THRESHOLD) {
            MainResult = 1;
        } else if (ANNResults[0] < CARC_THRESHOLD) {
            MainResult = 0;
        } else {
            MainResult = -1;
        }
        CurOutput.setMainResultValue(MainResult);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // Carcinogenicity classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for carcinogenicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }
        Res[1] = Format_3D.format(ANNResults[0]); // Prob for carcinogen
        Res[2] = Format_3D.format(ANNResults[1]); // Prob for non-carcinogen

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
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

        // Sets NN reliability index
        ADIndexNNReliability adRel = new ADIndexNNReliability();
        adRel.SetDelta(Math.abs(ANNResults[0] - ANNResults[1]));
        CurOutput.addADIndex(adRel);

        // Sets Neuron concordance index
        int nX = (int)ANNResults[2];
        int nY = (int)ANNResults[3];
        int pred = (int)CurOutput.getMainResultValue();
        NeuronChecker.Check(nX, nY, pred);
        ADIndexNNNeuronsConcordance adNNConc = new ADIndexNNNeuronsConcordance();
        adNNConc.SetNeuronData(NeuronChecker.Population, NeuronChecker.Concordance);
        CurOutput.addADIndex(adNNConc);

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
        double nnContribution = adNNConc.GetIndexValue();
        double deltaContribution = adRel.GetAssessmentClass()== ADIndex.INDEX_LOW ? 0 : 1;
        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution *
                nnContribution * deltaContribution;

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
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);

    }

}
