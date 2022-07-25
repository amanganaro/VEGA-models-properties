package insilico.sqfu;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.sqfu.descriptors.EmbeddedDescriptors;
import lombok.extern.log4j.Log4j;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

@Log4j
public class ismSQfu extends InsilicoModel {

    private static final String ModelData = "/data/model_sqfu.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;

    public ismSQfu() throws InitFailureException {

        super(ModelData);
        // Init PMML model
        try {
            URL src = getClass().getResource("/data/SQfu tree ensemble.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 20;
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
        this.DescriptorsNames[9] = "CATS2D_01_LL";
        this.DescriptorsNames[10] = "nCar";
        this.DescriptorsNames[11] = "SpMin1_Bh.i.";
        this.DescriptorsNames[12] = "Eta_betaP_A";
        this.DescriptorsNames[13] = "SM12_AEA.ri.";
        this.DescriptorsNames[14] = "N.";
        this.DescriptorsNames[15] = "nN.";
        this.DescriptorsNames[16] = "totalcharge";
        this.DescriptorsNames[17] = "CATS2D_00_PP";
        this.DescriptorsNames[18] = "SpMax2_Bh.m.";
        this.DescriptorsNames[19] = "C.024";



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

            Descriptors[0] = embeddedDescriptors.Alogp;
            Descriptors[1] = embeddedDescriptors.P_VSA_i_2;
            Descriptors[2] = embeddedDescriptors.MLogp;
            Descriptors[3] = embeddedDescriptors.P_VSA_p_3;
            Descriptors[4] = embeddedDescriptors.Eta_betaP;
            Descriptors[5] = embeddedDescriptors.CATS2D_00_LL;
            Descriptors[6] = embeddedDescriptors.C;
            Descriptors[7] = embeddedDescriptors.PCD;
            Descriptors[8] = embeddedDescriptors.Ui;
            Descriptors[9] = embeddedDescriptors.CATS2D_01_LL;
            Descriptors[10] = embeddedDescriptors.nCar;
            Descriptors[11] = embeddedDescriptors.SpMin1_Bh_i;
            Descriptors[12] = embeddedDescriptors.Eta_betaP_A;
            Descriptors[13] = embeddedDescriptors.SM12_AEA_ri;
            Descriptors[14] = embeddedDescriptors.N;
            Descriptors[15] = embeddedDescriptors.nN;
            Descriptors[16] = embeddedDescriptors.totalcharge;
            Descriptors[17] = embeddedDescriptors.CATS2D_00_PP;
            Descriptors[18] = embeddedDescriptors.SpMax2_Bh_m;
            Descriptors[19] = embeddedDescriptors.C_024;

            this.ExperimentalValue = embeddedDescriptors.ExperimentalValue;

        } catch (Throwable e) {
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
//        double ConvertedValue = Math.pow( (Prediction * 0.03 + 1), (1.0 / 0.03) ) * MW;
//        if (ConvertedValue>1)
//            Res[1] = Format_2D.format(ConvertedValue); // mg/L
//        else
//            Res[1] = Format_4D.format(ConvertedValue); // mg/L
//        Res[2] = Format_2D.format(MW); // MW
//        Res[3] = "-";
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
