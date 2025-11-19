package insilico.steroidogenesisedscreen.model.modelRFpmml;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.InitFailureException;
import org.dmg.pmml.FieldName;
import org.jpmml.evaluator.Evaluator;
import org.jpmml.evaluator.LoadingModelEvaluatorBuilder;
import org.jpmml.evaluator.mining.MiningProbabilityDistribution;

import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;


public class SteroidogenesisRFmodel {

    private final static String PMML_PATH = "/data/stero_eds_rf_model.pmml";
    private final SteroDescriptors descriptors;
    private final Evaluator evaluator;


    public SteroidogenesisRFmodel() throws InitFailureException {

        // init serialized model
        try {
            URL DataSource = SteroidogenesisRFmodel.class.getResource(PMML_PATH);
            evaluator = new LoadingModelEvaluatorBuilder()
                    .load(DataSource.openStream())
                    .build();
        } catch (Exception e) {
            throw new InitFailureException("Unable to load RF serialized model - " + e.getMessage());
        }

        descriptors = new SteroDescriptors();
    }


    public RFResult predict(String SMILES) {

        RFResult res = new RFResult();

        double[] features;
        try {
            features = descriptors.Scale(descriptors.Calculate(SMILES));
        } catch (Exception e) {
            return res;
        }

        Map<FieldName, Object> argumentsObject = new LinkedHashMap<>();
        for (int i=0; i<features.length; i++)
            argumentsObject.put(FieldName.create(SteroDescriptors.DESCRIPTOR_NAMES[i]), features[i]);

        try {
            Map<FieldName, ?> results = evaluator.evaluate(argumentsObject);

            MiningProbabilityDistribution mpd = (MiningProbabilityDistribution) results.get(FieldName.create("Target"));

            String mainPred = (String) mpd.getPrediction();
            if (mainPred.equalsIgnoreCase("0")) {
                res.Prediction = 0;
                if ( mpd.getProbability("0") != null )
                    res.Confidence = mpd.getProbability("0");
            } else if (mainPred.equalsIgnoreCase("1")) {
                res.Prediction = 1;
                if (mpd.getProbability("1") != null)
                    res.Confidence = mpd.getProbability("1");
            } else {
                res.Prediction = Descriptor.MISSING_VALUE;
                res.Confidence = 0;
            }

        } catch (Exception e) {
            return res;
        }

        return res;
    }

}
