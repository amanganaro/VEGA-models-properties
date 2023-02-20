package insilico.vapour_pressure;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.vapour_pressure.descriptors.EmbeddedDescriptors;
import insilico.vapour_pressure.runner.nrNetwork;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismVapourPressure extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismVapourPressure.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_vapour_pressure.xml";

    public ismVapourPressure() throws InitFailureException {
        super(ModelData);
        this.DescriptorsSize = 37;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "H-050";
        DescriptorsNames[1] = "F02[F-F]";
        DescriptorsNames[2] = "SM02_AEA(dm)";
        DescriptorsNames[3] = "SpMax_AEA(ed)";
        DescriptorsNames[4] = "ATS2m";
        DescriptorsNames[5] = "P_VSA_PPP_A";
        DescriptorsNames[6] = "NPerc";
        DescriptorsNames[7] = "piID";
        DescriptorsNames[8] = "CATS2D_07_DD";
        DescriptorsNames[9] = "CATS2D_01_DN";
        DescriptorsNames[10] = "N-072";
        DescriptorsNames[11] = "nAmideE-";
        DescriptorsNames[12] = "X3v";
        DescriptorsNames[13] = "TI1_L";
        DescriptorsNames[14] = "SpMaxA_Dz(m)";
        DescriptorsNames[15] = "BB_SA12";
        DescriptorsNames[16] = "nCONN";
        DescriptorsNames[17] = "CATS2D_07_DL";
        DescriptorsNames[18] = "B03[F-F]";
        DescriptorsNames[19] = "SaaNH";
        DescriptorsNames[20] = "ED_3";
        DescriptorsNames[21] = "CATS2D_09_DA";
        DescriptorsNames[22] = "LogP";
        DescriptorsNames[23] = "nRNH2";
        DescriptorsNames[24] = "F07[C-S]";
        DescriptorsNames[25] = "F10[C-O]";
        DescriptorsNames[26] = "CATS2D_09_AL";
        DescriptorsNames[27] = "BB_SA31c";
        DescriptorsNames[28] = "SpMin5_Bh(m)";
        DescriptorsNames[29] = "X1Av";
        DescriptorsNames[30] = "N-068";
        DescriptorsNames[31] = "nRNR2";
        DescriptorsNames[32] = "F09[O-O]";
        DescriptorsNames[33] = "IC1";
        DescriptorsNames[34] = "B01[C-N]";
        DescriptorsNames[35] = "nN=N";
        DescriptorsNames[36] = "F10[C-Cl]";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Vapour Pressure [log10(atm)]";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors EmbeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = EmbeddedDescriptors.H_050;
            Descriptors[1] = EmbeddedDescriptors.F02_F_F;
            Descriptors[2] = EmbeddedDescriptors.SM02_AEAdm;
            Descriptors[3] = EmbeddedDescriptors.SpMax_AEAed;
            Descriptors[4] = EmbeddedDescriptors.ATS2m;
            Descriptors[5] = EmbeddedDescriptors.P_VSA_PPP_A;
            Descriptors[6] = EmbeddedDescriptors.NPerc;
            Descriptors[7] = EmbeddedDescriptors.piID;
            Descriptors[8] = EmbeddedDescriptors.CATS2D_07_DD;
            Descriptors[9] = EmbeddedDescriptors.CATS2D_01_DN;
            Descriptors[10] = EmbeddedDescriptors.N_072;
            Descriptors[11] = EmbeddedDescriptors.nAmideE_;
            Descriptors[12] = EmbeddedDescriptors.X3v;
            Descriptors[13] = EmbeddedDescriptors.TI1_L;
            Descriptors[14] = EmbeddedDescriptors.SpMaxA_Dzm;
            Descriptors[15] = EmbeddedDescriptors.BB_SA12;
            Descriptors[16] = EmbeddedDescriptors.nCONN;
            Descriptors[17] = EmbeddedDescriptors.CATS2D_07_DL;
            Descriptors[18] = EmbeddedDescriptors.B03_F_F;
            Descriptors[19] = EmbeddedDescriptors.SaaNH;
            Descriptors[20] = EmbeddedDescriptors.ED_3;
            Descriptors[21] = EmbeddedDescriptors.CATS2D_09_DA;
            Descriptors[22] = EmbeddedDescriptors.LogP;
            Descriptors[23] = EmbeddedDescriptors.nRNH2;
            Descriptors[24] = EmbeddedDescriptors.F07_C_S;
            Descriptors[25] = EmbeddedDescriptors.F10_C_O;
            Descriptors[26] = EmbeddedDescriptors.CATS2D_09_AL;
            Descriptors[27] = EmbeddedDescriptors.BB_SA31c;
            Descriptors[28] = EmbeddedDescriptors.SpMin5_Bhm;
            Descriptors[29] = EmbeddedDescriptors.X1Av;
            Descriptors[30] = EmbeddedDescriptors.N_068;
            Descriptors[31] = EmbeddedDescriptors.nRNR2;
            Descriptors[32] = EmbeddedDescriptors.F09_O_O;
            Descriptors[33] = EmbeddedDescriptors.IC1;
            Descriptors[34] = EmbeddedDescriptors.B01_C_N;
            Descriptors[35] = EmbeddedDescriptors.nN_N;
            Descriptors[36] = EmbeddedDescriptors.F10_C_Cl;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        double Prediction;
        try {
            nrNetwork nn = nrNetwork.ReadFromFile("VegaModels-VapourPressure\\src\\main\\resources\\vapourPressure.nn");
            Prediction = nn.Calculate(Descriptors, true);
            CurOutput.setMainResultValue(Prediction);

            for (int i=0; i<Descriptors.length; i++){
                System.out.print("\t" + Descriptors[i]);
            }
            System.out.println();

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
