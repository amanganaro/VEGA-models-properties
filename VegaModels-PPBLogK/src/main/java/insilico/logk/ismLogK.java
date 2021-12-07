package insilico.logk;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.logk.descriptors.EmbeddedDescriptors;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

@Slf4j
public class ismLogK extends InsilicoModel {

    private static final String ModelData = "/data/model_logk.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;


    public ismLogK() throws InitFailureException {
        super(ModelData);

        try {
            URL src = getClass().getResource("/data/logk tree ensemble.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 24;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "ALOGP";
        this.DescriptorsNames[1] = "P_VSA_i_2";
        this.DescriptorsNames[2] = "MLOGP";
        this.DescriptorsNames[3] = "P_VSA_p_3";
        this.DescriptorsNames[4] = "Eta_betaP";
        this.DescriptorsNames[5] = "CATS2D_00_LL";
        this.DescriptorsNames[6] = "C.";
        this.DescriptorsNames[7] = "PCD";
        this.DescriptorsNames[8] = "Ui";
        this.DescriptorsNames[9] = "N.";
        this.DescriptorsNames[10] = "totalcharge";
        this.DescriptorsNames[11] = "CATS2D_00_PP";
        this.DescriptorsNames[12] = "C.024";
        this.DescriptorsNames[13] = "nCsp2";
        this.DescriptorsNames[14] = "MLOGP2";
        this.DescriptorsNames[15] = "GATS1i";
        this.DescriptorsNames[16] = "SpMax2_Bh.p.";
        this.DescriptorsNames[17] = "nBM";
        this.DescriptorsNames[18] = "MATS5e";
        this.DescriptorsNames[19] = "AMW";
        this.DescriptorsNames[20] = "F01.C.N.";
        this.DescriptorsNames[21] = "T.O..O.";
        this.DescriptorsNames[22] = "J_D.Dt";
        this.DescriptorsNames[23] = "SpMax_AEA.dm.";

        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Prediction";
        this.ResultsName[1] = "Experimental value [mg/l]";
    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule, false);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.ALogP;
            Descriptors[1] = embeddedDescriptors.P_VSA_i_2;
            Descriptors[2] = embeddedDescriptors.MLogP;
            Descriptors[3] = embeddedDescriptors.P_VSA_p_3;
            Descriptors[4] = embeddedDescriptors.Eta_betaP;
            Descriptors[5] = embeddedDescriptors.CATS2D_00_LL;
            Descriptors[6] = embeddedDescriptors.C;
            Descriptors[7] = embeddedDescriptors.PCD;
            Descriptors[8] = embeddedDescriptors.Ui;
            Descriptors[9] = embeddedDescriptors.N;
            Descriptors[10] = embeddedDescriptors.totalcharge;
            Descriptors[11] = embeddedDescriptors.CATS2D_00_PP;
            Descriptors[12] = embeddedDescriptors.C_024;
            Descriptors[13] = embeddedDescriptors.nCsp2;
            Descriptors[14] = embeddedDescriptors.MLOGP2;
            Descriptors[15] = embeddedDescriptors.GATS1i;
            Descriptors[16] = embeddedDescriptors.SpMax2_Bh_P;
            Descriptors[17] = embeddedDescriptors.nBm;
            Descriptors[18] = embeddedDescriptors.MATS5e;
            Descriptors[19] = embeddedDescriptors.AMW;
            Descriptors[20] = embeddedDescriptors.F01_C_N;
            Descriptors[21] = embeddedDescriptors.T_O_O;
            Descriptors[22] = embeddedDescriptors.J_D_DT;
            Descriptors[23] = embeddedDescriptors.SpMax_AEA_dm;

            this.ExperimentalValue = embeddedDescriptors.ExperimentalValue;



        } catch (Throwable e){
            return DESCRIPTORS_ERROR;
        }
        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        Map<String, Object> argumentsObject = new LinkedHashMap<>();
        try {
            for (int i=0; i<DescriptorsSize; i++)
                argumentsObject.put(this.DescriptorsNames[i], Descriptors[i]);
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

        // Run pmml model
        double Prediction;
        try {
            Prediction = Model.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Prediction);
        Res[1] = String.valueOf(ExperimentalValue);

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

    }
}
