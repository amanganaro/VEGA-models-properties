package insilico.dio1.model.modelGBTpmml;

import insilico.core.descriptor.Descriptor;

public class GBTResult {
    public String SMILES = "";
    public int prediction = Descriptor.MISSING_VALUE;
    public double probability = 0;
    public double confidence = 0;

}
