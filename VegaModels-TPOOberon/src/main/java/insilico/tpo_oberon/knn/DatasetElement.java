package insilico.tpo_oberon.knn;

import insilico.core.molecule.InsilicoMolecule;

public class DatasetElement {

//    public static final short MOLECULE_TRAINING =1;
//    public static final short MOLECULE_TEST =2;

    private InsilicoMolecule inputMolecule;
    private double[] descriptors;
    private int experimentalValue;
    private int status;


    public DatasetElement(InsilicoMolecule inputMolecule, double[] descriptors, int experimentalValue, int status) {
        this.inputMolecule = inputMolecule;
        this.descriptors = descriptors;
        this.experimentalValue = experimentalValue;
        this.status = status;
    }

    public DatasetElement(InsilicoMolecule inputMolecule, double[] descriptors) {
        this.inputMolecule = inputMolecule;
        this.descriptors = descriptors;
        this.experimentalValue = -999;
        this.status = -1;
    }

    public InsilicoMolecule getInputMolecule() {
        return inputMolecule;
    }

    public void setInputMolecule(InsilicoMolecule inputMolecule) {
        this.inputMolecule = inputMolecule;
    }

    public double[] getDescriptors() {
        return descriptors;
    }

    public void setDescriptors(double[] descriptors) {
        this.descriptors = descriptors;
    }

    public int getExperimentalValue() {
        return experimentalValue;
    }

    public void setExperimentalValue(int experimentalValue) {
        this.experimentalValue = experimentalValue;
    }

    public int getStatus() {
        return status;
    }

    public void setStatus(short status) {
        this.status = status;
    }

}
