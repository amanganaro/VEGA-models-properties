package insilico.daphnia_demetra;

public abstract class Predictor {

    String[] descrName, descrDetails;
    String measureUnit;
    String secondMeasureUnit;

    double min = 0.0, max = 0.0;

    public double getMin(){
        return min;
    }

    public double getMax(){
        return max;
    }

    public String getSecondMeasureUnit(){
        return secondMeasureUnit;
    }

    public String getMeasureUnit(){
        return measureUnit;
    }

    public int getDescriptorsNumber(){
        return descrName.length;
    }

    public String[] getDescriptorsNamesList(){
        return descrName;
    }

    public String[] getDescriptorsDetailsList(){
        return descrDetails;
    }

    public abstract double predict(double[] input);
    public abstract double minimumSubmodelsPrediction(double[] input);
    public abstract double maximumSubmodelsPrediction(double[] input);


}
