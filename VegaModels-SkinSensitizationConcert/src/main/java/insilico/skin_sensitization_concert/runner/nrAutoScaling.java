package insilico.skin_sensitization_concert.runner;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import insilico.core.descriptor.Descriptor;

import java.io.Serializable;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@JsonAutoDetect(fieldVisibility = JsonAutoDetect.Visibility.ANY, getterVisibility = JsonAutoDetect.Visibility.NONE, setterVisibility = JsonAutoDetect.Visibility.NONE)
public class nrAutoScaling implements Serializable {

    private double Mean;
    private double Variance;


    public nrAutoScaling(double Mean, double Variance) {
        this.Mean = Mean;
        this.Variance = Variance;
    }

    public nrAutoScaling() {
        this(0,1);
    }

    public double Scale(double Value) {

        if (Double.isNaN(Value))
            return Value;
        if (Double.isInfinite(Value))
            return Value;
        if (Value == Descriptor.MISSING_VALUE)
            return Value;

        // scaled values are rounded to 4 decimal digits
        double val = (Value - Mean) / Variance;
        val = Math.round(val * 10000) / 10000.0;
        return val;
    }

    public double Unscale(double Value) {

        if (Double.isNaN(Value))
            return Value;
        if (Double.isInfinite(Value))
            return Value;
        if (Value == Descriptor.MISSING_VALUE)
            return Value;

        return (Value * Variance) + Mean;
    }

    public nrAutoScaling getCopy() {
        return new nrAutoScaling(this.Mean, this.Variance);
    }

    /**
     * @return the Mean
     */
    public double getMean() {
        return Mean;
    }

    /**
     * @return the Variance
     */
    public double getVariance() {
        return Variance;
    }

    public void setMean(double mean) {
        Mean = mean;
    }

    public void setVariance(double variance) {
        Variance = variance;
    }
}
