package insilico.carcinogenicity_rat_female;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.coral.CoralModel;
import insilico.core.coral.models.carcino.CoralFRCarcinogenicity;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;

@Slf4j
public class ismCarcinogenicityRatFemale extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_carcinogenicity_female.xml";

    private CoralModel CarcinogenicityFemaleRat;

    public ismCarcinogenicityRatFemale() throws InitFailureException {
        super(ModelData);

        try {
            CarcinogenicityFemaleRat = new CoralFRCarcinogenicity();
        } catch (Exception ex){
            throw new InitFailureException("Unable to init coral model");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Carcinogenicity in female rat";

        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();
    }

    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {

        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        return blocks;
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
        double Prediction;
        try {
            Prediction = CarcinogenicityFemaleRat.Predict(this.CurMolecule.getInputSMILES());
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        this.ResultsSize = 1;
        String[] Res = new String[ResultsSize];
        Res[0] = Format_4D.format(Prediction);


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

        // Add transformed (mg/L) experimental if needed
//        if (CurOutput.HasExperimental()) {
//            double ConvertedValue = Math.pow(10, CurOutput.getExperimental());
//            if (ConvertedValue>1)
//                CurOutput.getResults()[2] = Format_2D.format(ConvertedValue); // 1/(mg/kg-day)
//            else
//                CurOutput.getResults()[2] = Format_4D.format(ConvertedValue); // 1/(mg/kg-day)
//        }

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
