package insilico.skin_sensitization_concert;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.skin_sensitization_concert.descriptors.EmbeddedDescriptors;
import insilico.skin_sensitization_concert.runner.nrNetwork;

import java.io.DataInputStream;
import java.net.URL;

public class ismSkinSensitizationConcert extends InsilicoModel {


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization_concert.xml";
    private static final String NNData = "/data/skinSensitization.nn";

    public ismSkinSensitizationConcert() throws InitFailureException {
        super(ModelData);
        this.DescriptorsSize = 17;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "Eig12_AEA(ri)";
        DescriptorsNames[1] = "nArOR";
        DescriptorsNames[2] = "nRCONH2";
        DescriptorsNames[3] = "nROH";
        DescriptorsNames[4] = "CATS2D_08_DA";
        DescriptorsNames[5] = "MATS3e";
        DescriptorsNames[6] = "N-075";
        DescriptorsNames[7] = "F05[C-O]";
        DescriptorsNames[8] = "H-054";
        DescriptorsNames[9] = "nR=Cs";
        DescriptorsNames[10] = "P_VSA_LogP_8";
        DescriptorsNames[11] = "CATS2D_09_AP";
        DescriptorsNames[12] = "SRW03";
        DescriptorsNames[13] = "nSH";
        DescriptorsNames[14] = "F06[N-F]";
        DescriptorsNames[15] = "B06[C-O]";
        DescriptorsNames[16] = "BB_SA18";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted skin sensitization";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors EmbeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = EmbeddedDescriptors.Eig12_AEA_ri;
            Descriptors[1] = EmbeddedDescriptors.nArOR;
            Descriptors[2] = EmbeddedDescriptors.nRCONH2;
            Descriptors[3] = EmbeddedDescriptors.nROH;
            Descriptors[4] = EmbeddedDescriptors.CATS2D_08_DA;
            Descriptors[5] = EmbeddedDescriptors.MATS3e;
            Descriptors[6] = EmbeddedDescriptors.N_075;
            Descriptors[7] = EmbeddedDescriptors.F05_C_O;
            Descriptors[8] = EmbeddedDescriptors.H_054;
            Descriptors[9] = EmbeddedDescriptors.nR_Cs;
            Descriptors[10] = EmbeddedDescriptors.P_VSA_LogP_8;
            Descriptors[11] = EmbeddedDescriptors.CATS2D_09_AP;
            Descriptors[12] = EmbeddedDescriptors.SRW03;
            Descriptors[13] = EmbeddedDescriptors.nSH;
            Descriptors[14] = EmbeddedDescriptors.F06_N_F;
            Descriptors[15] = EmbeddedDescriptors.B06_C_O;
            Descriptors[16] = EmbeddedDescriptors.BB_SA18;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        double Prediction;
        try {
            DataInputStream in;
            URL nnURL = getClass().getResource(NNData);
            in = new DataInputStream(nnURL.openStream());
            nrNetwork nn = nrNetwork.ReadFromFile(in);

            Prediction = nn.Calculate(Descriptors, true);
            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction);
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex) {
            return MODEL_ERROR;
        }
    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.AddMappingToNegativeValue(-1);
        adq.setMoleculesForIndexSize(3);
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

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.9, 0.65);
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
