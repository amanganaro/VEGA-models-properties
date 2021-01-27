package insilico.persistence_soil_irfmn;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAPersistenceSoil;
import insilico.core.constant.InsilicoConstants;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.knn.insilicoKnnQualitative;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;


import java.util.ArrayList;

import static insilico.core.constant.InsilicoConstants.KEY_ALERT_PERS_SOIL_NP;
import static insilico.core.constant.InsilicoConstants.KEY_ALERT_PERS_SOIL_VP;

/**
 *
 * @author User
 */
@Slf4j
public class ismPersistenceSoilIrfmn extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_pers_soil_irfmn.xml";

    private final String PersistenceAlertSet;
    private final insilicoKnnQualitative KNN;
    private insilicoKnnPrediction KnnPrediction;


    public ismPersistenceSoilIrfmn()
            throws InitFailureException {
        super(ModelData);

        // Build SA list
        PersistenceAlertSet = AlertEncoding.MergeAlertIds(new SAPersistenceSoil().getAlerts());

        // Build model object
        KNN =  new insilicoKnnQualitative();
        KNN.setNeighboursNumber(3);
        KNN.setMinSimilarity(0.7);
        KNN.setMinSimilarityForSingleResult(0.9);
        KNN.setEnhanceWeightFactor(2);
        KNN.setUseExperimentalRange(false);

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted persistence in soil";
        this.ResultsName[1] = "Molecules used for prediction";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexAlertsPersistenceWater().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();

    }


    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
        blocks.add(InsilicoConstants.SA_BLOCK_PERSISTENCE_SOIL_IRFMN);
        return blocks;
    }


    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(PersistenceAlertSet, a.getId())) {
                FoundSAs.add((Alert)a.clone());
            }
        return FoundSAs;
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

        // if molecule is not predicted, try to perform prediction with alerts
        if (KnnPrediction.getStatus() == insilicoKnnPrediction.KNN_MISSING_NO_MOLECULES) {
            try {
                AlertList FoundSAs = GetCalculatedAlert();
                if (FoundSAs.size() > 0) {
                    int nNP=0, nVP=0;
                    for (Alert CurAlert : FoundSAs.getSAList()) {
                        if (CurAlert.getBoolProperty(KEY_ALERT_PERS_SOIL_VP))
                            nVP++;
                        if (CurAlert.getBoolProperty(KEY_ALERT_PERS_SOIL_NP))
                            nNP++;
                    }
                    if ((nVP > 0)&&(nNP == 0))
                        KnnPrediction.setPrediction(4);
                    if ((nVP == 0)&&(nNP > 0))
                        KnnPrediction.setPrediction(1);
                }
            } catch (CloneNotSupportedException ex) {
                // do nothing,just skip SA check
            }
        }

        CurOutput.setMainResultValue(KnnPrediction.getPrediction());

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(KnnPrediction.getPrediction()); // Persistence classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for persistence value " + KnnPrediction.getPrediction());
            Res[0] = Double.toString(KnnPrediction.getPrediction());
        }
        Res[1] = String.valueOf(KnnPrediction.getNeighbours().size());
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        // Ad is performed on the K molecules used for the KNN model
        // or on 3 molecules if model used no molecules (only alerts)
        int UsedMols = KnnPrediction.getNeighbours().isEmpty() ?
                3 : KnnPrediction.getNeighbours().size();
        adq.setMoleculesForIndexSize(UsedMols);

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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.9, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.9, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets alerts contribution
        try {
            ADIndexAlertsPersistenceSoil adalerts = new ADIndexAlertsPersistenceSoil();
            adalerts.SetAlerts(GetCalculatedAlert(), CurOutput, KnnPrediction.getNeighbours().size());
            CurOutput.addADIndex(adalerts);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double saContribution = CurOutput.getADIndex(ADIndexAlertsPersistenceSoil.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution * saContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.85, 0.65);
        CurOutput.setADI(ADI);

        // Sets SA
        try {
            ADCheckSA ADSA = new ADCheckSA(TS);
            if (!ADSA.Calculate(CurMolecule, GetCalculatedAlert(), CurOutput, adq.GetSimilarMolecules()))
                return InsilicoModel.AD_ERROR;
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        if (CurOutput.getMainResultValue() == Descriptor.MISSING_VALUE) {

            CurOutput.setAssessment("N/A");
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_NA, "N/A"));

        } else {

            ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        }

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == Descriptor.MISSING_VALUE)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else {
            if (Val == 1)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
            else if (Val == 2)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
            else if (Val == 3)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
            else if (Val == 4)
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
            else
                CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
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
