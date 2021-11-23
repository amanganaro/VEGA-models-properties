package insilico.pparg_up;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.pparg_up.descriptors.EmbeddedDescriptors;
import lombok.extern.slf4j.Slf4j;
import org.dmg.pmml.FieldName;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

@Slf4j
public class ismPPARGup extends InsilicoModel {

    private static final String ModelData = "/data/model_ppargup.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;

    private static final double[] mean = {1.4872400881,0.8144922907,97.0816365639,59.0704446061,8.2034140969,0.3744493392,0.904185022,2.2676211454,8.1162301762,2.3448491189,1.1294581498,0.9471655026,1.1183964758,0.8032577093,3.9669603524,101.9245792952,2.8039647577,7.1685022026,38.7528634361,1.204845815,0.6057268722,18.3764845815};
    private static final double[] stdDeviation = {0.1558724554,0.1853956791,51.024656699,44.5739002958,6.1125838911,0.7488964499,1.4003717009,2.1057929397,8.6603366519,1.6189180219,0.0179934268,0.9851942441,0.2856590462,0.2013561017,5.6166598158,126.4102967879,3.8458594305,8.3449090676,6.9572107347,1.4203974312,0.48896333,25.4342756442};

    public ismPPARGup() throws InitFailureException {
        super(ModelData);

        try {
            URL src = getClass().getResource("/data/PPARg_up_RF_model.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Predicted Exp");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 22;
        this.DescriptorsNames = new String[DescriptorsSize];

        this.DescriptorsNames[0] = "SpMAD_B.m.";
        this.DescriptorsNames[1] = "rGes";
        this.DescriptorsNames[2] = "P_VSA_v_3";
        this.DescriptorsNames[3] = "P_VSA_LogP_4";
        this.DescriptorsNames[4] = "O.";
        this.DescriptorsNames[5] = "O.057";
        this.DescriptorsNames[6] = "nCconj";
        this.DescriptorsNames[7] = "nCb.";
        this.DescriptorsNames[8] = "MLOGP2";
        this.DescriptorsNames[9] = "MLOGP";
        this.DescriptorsNames[10] = "Mi";
        this.DescriptorsNames[11] = "GATS8s";
        this.DescriptorsNames[12] = "GATS1p";
        this.DescriptorsNames[13] = "GATS1m";
        this.DescriptorsNames[14] = "F06.C.O.";
        this.DescriptorsNames[15] = "D.Dtr06";
        this.DescriptorsNames[16] = "CATS2D_07_AL";
        this.DescriptorsNames[17] = "CATS2D_03_LL";
        this.DescriptorsNames[18] = "C.";
        this.DescriptorsNames[19] = "C.026";
        this.DescriptorsNames[20] = "B07.C.O.";
        this.DescriptorsNames[21] = "ATSC7m";

        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Prediction";
        this.ResultsName[1] = "0 Probability";
        this.ResultsName[2] = "1 Probability";
        this.ResultsName[3] = "Experimental value";

    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule, false);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.SpMAD_B_m;
            Descriptors[1] = embeddedDescriptors.rGes;
            Descriptors[2] = embeddedDescriptors.P_VSA_v_3;
            Descriptors[3] = embeddedDescriptors.P_VSA_LogP_4;
            Descriptors[4] = embeddedDescriptors.O;
            Descriptors[5] = embeddedDescriptors.O_057;
            Descriptors[6] = embeddedDescriptors.nCconj;
            Descriptors[7] = embeddedDescriptors.nCb;
            Descriptors[8] = embeddedDescriptors.MLOGP2;
            Descriptors[9] = embeddedDescriptors.MLOGP;
            Descriptors[10] = embeddedDescriptors.Mi;
            Descriptors[11] = embeddedDescriptors.GATS8s;
            Descriptors[12] = embeddedDescriptors.GATS1p;
            Descriptors[13] = embeddedDescriptors.GATS1m;
            Descriptors[14] = embeddedDescriptors.F06_C_O;
            Descriptors[15] = embeddedDescriptors.D_Dtr06;
            Descriptors[16] = embeddedDescriptors.CATS2D_07_AL;
            Descriptors[17] = embeddedDescriptors.CATS2D_03_LL;
            Descriptors[18] = embeddedDescriptors.C;
            Descriptors[19] = embeddedDescriptors.C_026;
            Descriptors[20] = embeddedDescriptors.B07_C_O;
            Descriptors[21] = embeddedDescriptors.ATSC7m;
            System.out.println();

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }



    @Override
    protected short CalculateModel() {
        Map<String, Object> argumentsObject = new LinkedHashMap<>();
        double[] ScaledDescriptors = new double[this.DescriptorsSize];

        for (int i=0; i<DescriptorsSize; i++) {
//            if(i == 3)
//                ScaledDescriptors[i] = Descriptors[i];
            ScaledDescriptors[i] = (Descriptors[i] - mean[i]) / stdDeviation[i];
        }
        try {
            for (int i=0; i<DescriptorsSize; i++)
                argumentsObject.put(this.DescriptorsNames[i], ScaledDescriptors[i]);
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


        // Run pmml model
        int Prediction;
        double Probability0;
        double Probability1;
        try {
//            Prediction = Model.Evaluate(argumentsObject);
            Map<FieldName, ?> outputs = Model.EvaluateFullOutput(argumentsObject);
            Prediction = Integer.parseInt((String) outputs.get(FieldName.create("Predicted_Exp")));
            Probability0 = (Double) outputs.get(FieldName.create("Probability_0"));
            Probability1 = (Double) outputs.get(FieldName.create("Probability_1"));

        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Prediction);
        Res[1] = String.valueOf(Probability0);
        Res[2] = String.valueOf(Probability1);
        Res[3] = String.valueOf(ExperimentalValue);
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
    protected void CalculateAssessment() {}
}
