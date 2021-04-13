package insilico.pgp_nic.descriptors;

import insilico.core.exception.InitFailureException;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;

/**
 * Class for the Counter-Propagation NN model
 */
public class modelCPANN {

    private final static String URL_WEIGHT = "/data/nn_weights.txt";
    private final static String URL_NORM = "/data/nn_normalization.txt";
    
    private CPANNNeuron[][] NN;
    private int NN_dimX, NN_dimY;
    private double[] Desc_Mean;
    private double[] Desc_Std;
    
    
    /**
     * 
     * @throws InitFailureException
     */
    public modelCPANN() throws InitFailureException  {
        
        try {

            // Loads rules for the model
            LoadRules();
            
            // Loads params for the normalization of descriptors
            LoadNormalizationParameters();

            // Loads training set
            // TO DO
            
        } catch (Exception e) {
            throw new InitFailureException();
        }
        
    }
    
    
    /**
     * 
     * @throws Exception
     */
    private void LoadRules() throws Exception {

	int n=0, j=0, k=0;
        int Idx=0;
        int Idx_X=0, Idx_Y=0;
        
        URL u;

        u = getClass().getResource(URL_WEIGHT);
        InputStream is = u.openStream();
        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String strLine;

        // First line: number of neurons and variables
        strLine = br.readLine();
        String[] parsedStr = strLine.split("\t");
        NN_dimX = Integer.valueOf(parsedStr[0]).intValue();
        NN_dimY = Integer.valueOf(parsedStr[1]).intValue();

        NN = new CPANNNeuron[NN_dimX][NN_dimY];
        Idx_X = 0;
        Idx_Y = 0;

        while ((strLine = br.readLine()) != null)   {

            // Gets number of tokens
            parsedStr = strLine.split("\t");

            // Adds the new line read
            double[] BufNN = new double[CPANNNeuron.DIM_P];
            double[] BufOut = new double[CPANNNeuron.DIM_OUT];
            
            for (int i=0; i<BufNN.length; i++)
                BufNN[i] = Double.valueOf(parsedStr[0 + i]).doubleValue();
            for (int i=0; i<BufOut.length; i++)
                BufOut[i] = Double.valueOf(parsedStr[BufNN.length + i]).doubleValue();
            
            NN[Idx_X][Idx_Y] = new CPANNNeuron(BufNN, BufOut);

            Idx_Y++;
            if (Idx_Y>(NN_dimY-1)) {
                Idx_X++;
                Idx_Y = 0;
            }
        }
    }    


    /**
     * 
     * @throws Exception
     */
    private void LoadNormalizationParameters() throws Exception {

	int n=0, j=0, k=0;
        int Idx=0;
        int Idx_X=0, Idx_Y=0;
        
        URL u;

        u = getClass().getResource(URL_NORM);
        InputStream is = u.openStream();
        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String strLine;

        Desc_Mean = new double[CPANNNeuron.DIM_P];
        Desc_Std = new double[CPANNNeuron.DIM_P];

        // Skips first line (headers and text)
        strLine = br.readLine();

        // Second line: std
        strLine = br.readLine();
        String[] parsedStr = strLine.split("\t");
        for (int i=0; i<Desc_Std.length; i++)
            Desc_Std[i] = Double.valueOf(parsedStr[1 + i]).doubleValue();

        // Third line: mean
        strLine = br.readLine();
        parsedStr = strLine.split("\t");
        for (int i=0; i<Desc_Mean.length; i++)
            Desc_Mean[i] = Double.valueOf(parsedStr[1 + i]).doubleValue();

    }    

    
    /**
     * 
     * @param Descriptors
     * @return
     */
    public double[] CalculatePrediction(double[] Descriptors) {
        
        double[] DescValues = new double[Descriptors.length];
        System.arraycopy(Descriptors, 0, DescValues, 0, Descriptors.length);

        double CurMin=1000;
        double val;
        int MinX=0, MinY=0;
        double[] Res = new double[6];   // Response: Y1 + Y2 + Y3 +
                                        // NeuronX + NeuronY + Distance
        
        // Normalizes descriptors
        for (int i=0; i<DescValues.length; i++)
            DescValues[i] = (DescValues[i] - Desc_Mean[i]) / Desc_Std[i];
        
        // Search neuron with minimum distance
        for (int i=0; i<NN_dimX; i++)
            for (int j=0; j<NN_dimY; j++) {
                if ((val=NN[i][j].CalculateDistance(DescValues))<CurMin) {
                    CurMin = val;
                    MinX = i;
                    MinY = j;
                }
            }

        double ED = Math.sqrt(CurMin);
        
        // Builds response
        Res[0] = NN[MinX][MinY].Out_Layer[0]; // Y1
        Res[1] = NN[MinX][MinY].Out_Layer[1]; // Y2
        Res[2] = NN[MinX][MinY].Out_Layer[2]; // Y2
        Res[3] = MinX; // Chosen neuron X
        Res[4] = MinY; // Chosen neuron Y
        Res[5] = ED; // distance from the center of the neuron (neuron weights)
        
        
        return Res;
    }
    
}
