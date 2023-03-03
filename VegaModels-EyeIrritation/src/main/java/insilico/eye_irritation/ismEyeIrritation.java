package insilico.eye_irritation;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.eye_irritation.descriptors.EmbeddedDescriptors;
import insilico.eye_irritation.runner.nrNetwork;
import javassist.runtime.Desc;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;

public class ismEyeIrritation extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismEyeIrritation.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_eye_irritation.xml";

    public ismEyeIrritation() throws InitFailureException {
        super(ModelData);

        // Define no. of descriptors
        this.DescriptorsSize = 50;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "SpMaxA_AEA(dm)";
        DescriptorsNames[1] = "B04[N-O]";
        DescriptorsNames[2] = "MATS8s";
        DescriptorsNames[3] = "P_VSA_e_5";
        DescriptorsNames[4] = "CATS2D_09_DA";
        DescriptorsNames[5] = "N-067";
        DescriptorsNames[6] = "Eig07_AEA(dm)";
        DescriptorsNames[7] = "nOxiranes";
        DescriptorsNames[8] = "CATS2D_02_DN";
        DescriptorsNames[9] = "MATS3s";
        DescriptorsNames[10] = "nR=Ct";
        DescriptorsNames[11] = "nR#CH/X";
        DescriptorsNames[12] = "CATS2D_06_PP";
        DescriptorsNames[13] = "B05[N-S]";
        DescriptorsNames[14] = "N-068";
        DescriptorsNames[15] = "F06[N-F]";
        DescriptorsNames[16] = "B10[C-S]";
        DescriptorsNames[17] = "B08[N-S]";
        DescriptorsNames[18] = "ATS5i";
        DescriptorsNames[19] = "ATSC1m";
        DescriptorsNames[20] = "B06[C-N]";
        DescriptorsNames[21] = "nCH2RX";
        DescriptorsNames[22] = "CATS2D_03_AP";
        DescriptorsNames[23] = "CATS2D_04_DN";
        DescriptorsNames[24] = "B04[C-O]";
        DescriptorsNames[25] = "nSO3OH";
        DescriptorsNames[26] = "B01[N-S]";
        DescriptorsNames[27] = "GATS4p";
        DescriptorsNames[28] = "MATS7m";
        DescriptorsNames[29] = "B02[C-C]";
        DescriptorsNames[30] = "SpDiam_EA(dm)";
        DescriptorsNames[31] = "nRCOOR";
        DescriptorsNames[32] = "B02[O-S]";
        DescriptorsNames[33] = "SpMaxA_EA(dm)";
        DescriptorsNames[34] = "B04[C-S]";
        DescriptorsNames[35] = "B10[C-N]";
        DescriptorsNames[36] = "P_VSA_m_2";
        DescriptorsNames[37] = "Eig05_EA(dm)";
        DescriptorsNames[38] = "F04[N-Br]";
        DescriptorsNames[39] = "B07[C-S]";
        DescriptorsNames[40] = "CATS2D_08_AN";
        DescriptorsNames[41] = "Cl-088";
        DescriptorsNames[42] = "CATS2D_01_DN";
        DescriptorsNames[43] = "B03[O-S]";
        DescriptorsNames[44] = "GATS5m";
        DescriptorsNames[45] = "SM15_AEA(dm)";
        DescriptorsNames[46] = "SM03_AEA(dm)";
        DescriptorsNames[47] = "F05[N-N]";
        DescriptorsNames[48] = "B02[O-O]";
        DescriptorsNames[49] = "B01[C-Br]";


        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Eye Irritation class";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.SpMaxA_AEA_dm;
            Descriptors[1] = embeddedDescriptors.B04_N_O;
            Descriptors[2] = embeddedDescriptors.MATS8s;
            Descriptors[3] = embeddedDescriptors.P_VSA_e_5;
            Descriptors[4] = embeddedDescriptors.CATS2D_09_DA;
            Descriptors[5] = embeddedDescriptors.N_067;
            Descriptors[6] = embeddedDescriptors.Eig07_AEA_dm;
            Descriptors[7] = embeddedDescriptors.nOxiranes;
            Descriptors[8] = embeddedDescriptors.CATS2D_02_DN;
            Descriptors[9] = embeddedDescriptors.MATS3s;
            Descriptors[10] = embeddedDescriptors.nR_Ct;
            Descriptors[11] = embeddedDescriptors.nR_CH_X;
            Descriptors[12] = embeddedDescriptors.CATS2D_06_PP;
            Descriptors[13] = embeddedDescriptors.B05_N_S;
            Descriptors[14] = embeddedDescriptors.N_068;
            Descriptors[15] = embeddedDescriptors.F06_N_F;
            Descriptors[16] = embeddedDescriptors.B10_C_S;
            Descriptors[17] = embeddedDescriptors.B08_N_S;
            Descriptors[18] = embeddedDescriptors.ATS5i;
            Descriptors[19] = embeddedDescriptors.ATSC1m;
            Descriptors[20] = embeddedDescriptors.B06_C_N;
            Descriptors[21] = embeddedDescriptors.nCH2RX;
            Descriptors[22] = embeddedDescriptors.CATS2D_03_AP;
            Descriptors[23] = embeddedDescriptors.CATS2D_04_DN;
            Descriptors[24] = embeddedDescriptors.B04_C_O;
            Descriptors[25] = embeddedDescriptors.nSO3OH;
            Descriptors[26] = embeddedDescriptors.B01_N_S;
            Descriptors[27] = embeddedDescriptors.GATS4p;
            Descriptors[28] = embeddedDescriptors.MATS7m;
            Descriptors[29] = embeddedDescriptors.B02_C_C;
            Descriptors[30] = embeddedDescriptors.SpDiam_EA_dm;
            Descriptors[31] = embeddedDescriptors.nRCOOR;
            Descriptors[32] = embeddedDescriptors.B02_O_S;
            Descriptors[33] = embeddedDescriptors.SpMaxA_EA_dm;
            Descriptors[34] = embeddedDescriptors.B04_C_S;
            Descriptors[35] = embeddedDescriptors.B10_C_N;
            Descriptors[36] = embeddedDescriptors.P_VSA_m_2;
            Descriptors[37] = embeddedDescriptors.Eig05_EA_dm;
            Descriptors[38] = embeddedDescriptors.F04_N_Br;
            Descriptors[39] = embeddedDescriptors.B07_C_S;
            Descriptors[40] = embeddedDescriptors.CATS2D_08_AN;
            Descriptors[41] = embeddedDescriptors.Cl_088;
            Descriptors[42] = embeddedDescriptors.CATS2D_01_DN;
            Descriptors[43] = embeddedDescriptors.B03_O_S;
            Descriptors[44] = embeddedDescriptors.GATS5m;
            Descriptors[45] = embeddedDescriptors.SM15_AEA_dm;
            Descriptors[46] = embeddedDescriptors.SM03_AEA_dm;
            Descriptors[47] = embeddedDescriptors.F05_N_N;
            Descriptors[48] = embeddedDescriptors.B02_O_O;
            Descriptors[49] = embeddedDescriptors.B01_C_Br;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        double Prediction;
        try {
            nrNetwork nn = nrNetwork.ReadFromFile("VegaModels-EyeIrritation\\src\\main\\resources\\eyeIrritation.nn");
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
        // todo Insert Applicability Domain here
        return 0;
    }

    @Override
    protected void CalculateAssessment() {
        // todo Insert Model Assessment here
    }
}
