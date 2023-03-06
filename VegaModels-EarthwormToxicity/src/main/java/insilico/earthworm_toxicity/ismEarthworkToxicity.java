package insilico.earthworm_toxicity;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.earthworm_toxicity.descriptors.EmbeddedDescriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismEarthworkToxicity extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismEarthworkToxicity.class);


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_earthworm_toxicity.xml";

    public ismEarthworkToxicity() throws InitFailureException {
        super(ModelData);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "pNOEC [-log(mg/Kg)]";
        this.ResultsName[1] = "NOEC [mg/Kg]";

        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "B06[N−O]";
        DescriptorsNames[1] = "SaasC";
        DescriptorsNames[2] = "B06[N−F]";
        DescriptorsNames[3] = "GATS1m";
        DescriptorsNames[4] = "VE1signB(e)";
        DescriptorsNames[5] = "B03[N−Cl]";
        DescriptorsNames[6] = "B08[N−Cl]";
        DescriptorsNames[7] = "S−107";


    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.B06_N_O;
            Descriptors[1] = embeddedDescriptors.SaasC;
            Descriptors[2] = embeddedDescriptors.B06_N_F;
            Descriptors[3] = embeddedDescriptors.GATS1m;
            Descriptors[4] = embeddedDescriptors.VE1signB_e;
            Descriptors[5] = embeddedDescriptors.B03_N_Cl;
            Descriptors[6] = embeddedDescriptors.B08_N_Cl;
            Descriptors[7] = embeddedDescriptors.S_107;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        double Prediction;
        try {

            Prediction = 1.062
                        - 1.389 * Descriptors[0]
                        - 0.156 * Descriptors[1]
                        + 1.29 * Descriptors[2]
                        - 1.9 * Descriptors[3]
                        - 2.425 * Descriptors[4]
                        - 0.821 * Descriptors[5]
                        + 0.225 * Descriptors[6]
                        + 1.233 * Descriptors[7];

            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = String.valueOf(Format_3D.format(Prediction));
            Res[1] = String.valueOf(Format_3D.format(Math.pow(10,-Prediction)));
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex){
            return MODEL_ERROR;
        }

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
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.6);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.6);
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

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
