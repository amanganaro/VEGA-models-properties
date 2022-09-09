package insilico.glucocorticoid_receptor;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.glucocorticoid_receptor.alert.SAGlucocorticoidReceptor;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class ismGlucocorticoidReceptor extends InsilicoModel {
    private static final Logger log = LogManager.getLogger(ismGlucocorticoidReceptor.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_glucocorticoid_receptor.xml";

    private SAGlucocorticoidReceptor saGlucocorticoidReceptor;

    public ismGlucocorticoidReceptor() throws InitFailureException {
        super(ModelData);
        try {
            saGlucocorticoidReceptor = new SAGlucocorticoidReceptor();
        } catch (Exception e) {
            throw new InitFailureException("Unable to init smarts - " + e.getMessage());
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Receptor Activity";

        // Define AD items
        this.ADItemsName = new String[0];
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
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
            saGlucocorticoidReceptor.Match(CurMolecule);
        } catch (Exception ex){
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(saGlucocorticoidReceptor.getMatches());
        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(saGlucocorticoidReceptor.getMatches()); // aromatase classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for ED value " + saGlucocorticoidReceptor.getMatches());
            Res[0] = Integer.toString(saGlucocorticoidReceptor.getMatches());
        }

        CurOutput.setResults(Res);
        return MODEL_CALCULATED;

    }

    @Override
    protected short CalculateAD() {
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToPositiveValue(2);
        adq.AddMappingToPositiveValue(3);
        adq.AddMappingToNegativeValue(0);
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

        if (CurOutput.getMainResultValue() == -1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else if (CurOutput.getMainResultValue() == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
    }
}
