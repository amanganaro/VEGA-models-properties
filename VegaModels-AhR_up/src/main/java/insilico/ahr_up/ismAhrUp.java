package insilico.ahr_up;

import insilico.ahr_up.descriptor.EmbeddedDescriptors;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.weights.TopologicalDistances;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.pmml.ModelANNFromPMML;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dmg.pmml.FieldName;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

public class ismAhrUp extends InsilicoModel {

    private static final double[] mean = {22.9688277512, 1.6291110048, 0.4232038278, 0.0201980861, 4.7118889952, 0.3668679426, 41.4590545455};
    private static final double[] stdDeviation = {21.873482939, 0.2498676299, 0.1641064496, 0.0154338905, 1.530393509, 0.5306807074, 19.2944565571};

    private static final String ModelData = "/data/model_ahr_up.xml";
    private ModelANNFromPMML Model;
    private double ExperimentalValue = 0;

    private static final Logger log = LogManager.getLogger(ismAhrUp.class);


    public ismAhrUp() throws InitFailureException {
        super(ModelData);

        try {

            URL src = getClass().getResource("/data/RF_model_ahrup.pmml");
            Model = new ModelANNFromPMML(src.openStream(), "Predicted Exp");
        } catch (IOException | InitFailureException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }

        this.DescriptorsSize = 7;
        this.DescriptorsNames = new String[DescriptorsSize];

        this.DescriptorsNames[0] = "ATSC5m";
        this.DescriptorsNames[1] = "SpMin3_Bh.v.";
        this.DescriptorsNames[2] = "SpMaxA_B.s.";
        this.DescriptorsNames[3] = "ChiA_D";
        this.DescriptorsNames[4] = "ATS7s";
        this.DescriptorsNames[5] = "ATSC7e";
        this.DescriptorsNames[6] = "Si";

        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Prediction";
        this.ResultsName[1] = "Prediction_0";
        this.ResultsName[2] = "Prediction_1";
    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {

        try {
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule, false);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.ATSC5m;
            Descriptors[1] = embeddedDescriptors.SpMin3_Bh_v;
            Descriptors[2] = embeddedDescriptors.SpMaxA_B_s;
            Descriptors[3] = embeddedDescriptors.ChiA_D;
            Descriptors[4] = embeddedDescriptors.ATS7s;
            Descriptors[5] = embeddedDescriptors.ATSC7e;
            Descriptors[6] = embeddedDescriptors.Si;

        } catch (Throwable e){
            log.warn(e.getMessage());
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        Map<String, Object> argumentsObject = new LinkedHashMap<>();
        double[] ScaledDescriptors = new double[this.DescriptorsSize];

        for (int i=0; i<DescriptorsSize; i++) {
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
        double prediction_0;
        double prediction_1;

        try {
//            Prediction = Model.Evaluate(argumentsObject);
            Map<FieldName, ?> outputs = Model.EvaluateFullOutput(argumentsObject);
            Prediction = Integer.parseInt((String) outputs.get(FieldName.create("Predicted_Exp")));
            prediction_0 = (Double) outputs.get(FieldName.create("Probability_0"));
            prediction_1 = (Double) outputs.get(FieldName.create("Probability_1"));
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Prediction);
        Res[1] = String.valueOf(prediction_0);
        Res[2] = String.valueOf(prediction_1);

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
