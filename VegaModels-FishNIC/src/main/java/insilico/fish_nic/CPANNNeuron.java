package insilico.fish_nic;


/**
 * Class for NN neuron
 */
public class CPANNNeuron {
    
    public final static int DIM_P = 20;      // dimension of the weight
    public final static int DIM_OUT = 1;    // dimension of the output layer weight
    
    public double[] NN_Layer;
    public double[] Out_Layer;
    
    
    public CPANNNeuron(double[] NN_Values, double[] Out_Values) {
        NN_Layer = NN_Values;
        Out_Layer = Out_Values;
    }
    
    
    public double CalculateDistance(double[] ObjValues) {
        
        double dist = 0;
        
        if (ObjValues.length!=DIM_P)
            return -999;
        
        // Distance calculated as sum of squares
        
        for (int i=0; i<DIM_P; i++) 
            dist += Math.pow(ObjValues[i] - NN_Layer[i],2);
        
        
        return dist;
    }

}
