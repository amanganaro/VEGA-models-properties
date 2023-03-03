package insilico.eye_irritation.runner;

import com.fasterxml.jackson.databind.ObjectMapper;
import insilico.core.descriptor.DescriptorBlock;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;

public class nrNetwork implements Serializable {

    public final static short ACT_LINEAR = 1;
    public final static short ACT_SIGMOID = 2;

    public String Description;
    public String[] DescriptorsName;

    public nrNeuron[] Neurons;
    public int OutputNeuron;
    public int[] InputNeurons;

    public nrAutoScaling[] InputScalers;

    public boolean Classification;

    public double Min;
    public double Max;

    private transient ArrayList<DescriptorBlock> DescBlocks;


    public nrNetwork() {

        DescriptorsName = new String[0];
        DescBlocks = new ArrayList<>();
    }

    private String PruneName(String desc) {

        int idx = desc.indexOf("_SQ");
        if (idx != -1)
            if (idx == desc.length()-3)
                return desc.substring(0,desc.length()-3);

        idx = desc.indexOf("_CU");
        if (idx != -1)
            if (idx == desc.length()-3)
                return desc.substring(0,desc.length()-3);

        return desc;
    }

    private boolean DescriptorIsPower(String desc, String pow) {
        int idx = desc.indexOf(pow);
        if (idx != -1)
            if (idx == desc.length()-3)
                return true;
        return false;
    }





    public double Calculate(double[] X) throws Exception {
        return Calculate(X, true);
    }

    public double Calculate(double[] X, boolean ScaleInput) throws Exception {

        if (X.length != (InputNeurons.length-1))
            throw new Exception("Wrong size of input");

        for (int i=0; i<X.length; i++) {
            if (ScaleInput)
                Neurons[InputNeurons[i]].setConstVal( InputScalers[i].Scale(X[i]) );
            else
                Neurons[InputNeurons[i]].setConstVal(X[i]);
        }

        double val = RecursiveCalculateNeuron(OutputNeuron);

        if (this.Classification) {
            // set binary class
            val = (val > 0.5) ? 1.0 : 0.0;
        } else {
            // unscale using range scaling
            val = (val * (Max - Min)) + Min;

        }
        return val;
    }

    private double RecursiveCalculateNeuron(int NeuronIdx) {

        nrNeuron n = Neurons[NeuronIdx];

        if (n.getWeights().length == 0)
            return n.getConstVal();

        double[] v = new double[n.getWeights().length];
        for (int i=0; i<n.getWeights().length; i++) {
            v[i] = RecursiveCalculateNeuron(n.getAntecedents()[i]);
        }

        double res = 0;

        for (int i=0; i<n.getWeights().length; i++)
            res += v[i] * n.getWeights()[i];

        double SigmoidCoeff = 1;
        if (n.getActivationType() == ACT_SIGMOID) {
            res = 1.0 / (1.0 + Math.exp(-1 * res * SigmoidCoeff));
        }
        return res;
    }


    public void SaveToFile(String FileName) throws IOException {
        ObjectMapper objectMapper = new ObjectMapper();
        objectMapper.writeValue(new File(FileName), this);
    }


    public static nrNetwork ReadFromFile(String FileName) throws IOException {
        ObjectMapper objectMapper = new ObjectMapper();
        return objectMapper.readValue(new File(FileName), nrNetwork.class);
    }

    public static nrNetwork ReadFromFile(InputStream source) throws IOException {
        ObjectMapper objectMapper = new ObjectMapper();
        return objectMapper.readValue(source, nrNetwork.class);
    }

}
