package insilico.skin_sensitization_concert.runner;

public class nrNeuron {

    public int[] Antecedents;
    public double[] Weights;
    public double ConstVal;
    public short ActivationType;

    public nrNeuron() {

    }


    public int[] getAntecedents() {
        return Antecedents;
    }

    public void setAntecedents(int[] antecedents) {
        Antecedents = antecedents;
    }

    public double[] getWeights() {
        return Weights;
    }

    public void setWeights(double[] weights) {
        Weights = weights;
    }

    public double getConstVal() {
        return ConstVal;
    }

    public void setConstVal(double constVal) {
        ConstVal = constVal;
    }

    public short getActivationType() {
        return ActivationType;
    }

    public void setActivationType(short activationType) {
        ActivationType = activationType;
    }
}
