package insilico.steroidogenesisedscreen.model;

import insilico.core.descriptor.Descriptor;

public class SteroidogenesisConsResult {

    public int Result = Descriptor.MISSING_VALUE;

    public int KNN_result = Descriptor.MISSING_VALUE;
    public double KNN_similarity = 0;
    public int KNN_n_positive = 0;
    public int KNN_n_negative = 0;

    public int RF_result = Descriptor.MISSING_VALUE;
    public double RF_confidence = 0;
}
