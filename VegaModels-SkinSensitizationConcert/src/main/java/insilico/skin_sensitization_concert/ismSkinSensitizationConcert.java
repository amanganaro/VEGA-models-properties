package insilico.skin_sensitization_concert;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.skin_sensitization_concert.descriptors.EmbeddedDescriptors;
import insilico.skin_sensitization_concert.runner.nrNetwork;

public class ismSkinSensitizationConcert extends InsilicoModel {


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization_concert.xml";

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
            nrNetwork nn = nrNetwork.ReadFromFile("VegaModels-SkinSensitizationConcert/src/main/resources/skinSensitization.nn");
            Prediction = nn.Calculate(Descriptors, true);
            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = String.valueOf(Format_3D.format(Prediction));
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex) {
            return MODEL_ERROR;
        }
    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

    }
}
