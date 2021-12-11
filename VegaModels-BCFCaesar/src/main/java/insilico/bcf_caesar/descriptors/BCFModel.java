package insilico.bcf_caesar.descriptors;

import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.exception.InitFailureException;

import java.net.URL;

/**
 *
 * @author User
 */
public class BCFModel {

    private Matrix GACentres;   // Centres for GA model
    private Matrix GAWeights;   // Weights for GA model
    private Matrix HMCentres;   // Centres for HM model
    private Matrix HMWeights;   // Weights for HM model
    private Matrix Ranges;      // Loaded ranges
    private Matrix GARanges;    // Ranges for GA descriptors
    private Matrix HMRanges;    // Ranges for HM descriptors

    private double Prediction;
    private double PredictionHM;
    private double PredictionGA;
    
    
    public BCFModel() throws InitFailureException {
        
        RBFNNDataReader DataCentresReader;
        RBFNNDataReader DataWeightsReader;
        RBFNNDataReader DataSubsetReader;
        RBFNNDataReader RangesReader;
        Matrix BufMatrix;
        int[][] DataSubset;

        try {
            
            URL u;
            
            u = getClass().getResource("/data/bcfcaesar_RBFNN_ranges.csv");
            RangesReader = new RBFNNDataReader(u.openStream(),0,0,0);
            Ranges = new Matrix(RangesReader.Data);

            // Reads data for GA model
            
            u = getClass().getResource("/data/bcfcaesar_GA_centres.csv");
            DataCentresReader = new RBFNNDataReader(u.openStream(),0,0,0);
            u = getClass().getResource("/data/bcfcaesar_GA_weights.csv");
            DataWeightsReader = new RBFNNDataReader(u.openStream(),0,0,0);
            u = getClass().getResource("/data/bcfcaesar_GA_subset.csv");
            DataSubsetReader = new RBFNNDataReader(u.openStream(),0,0,0);

            GACentres = new Matrix(DataCentresReader.Data);
            GAWeights = new Matrix(DataWeightsReader.Data);
            DataSubset = DataSubsetReader.DataAsInt();

            BufMatrix = new Matrix(GACentres.getRowDimension(), DataSubsetReader.dim_j);
            for (int i=0; i<GACentres.getRowDimension(); i++)
                for (int j=0; j<DataSubsetReader.dim_j; j++)
                    BufMatrix.set(i, j, GACentres.get(i, DataSubset[0][j]));
            GACentres = BufMatrix;

            GARanges = new Matrix(2,5);
            for (int j=0; j<5; j++) {
                for (int i=0; i<2; i++) 
                    GARanges.set(i, j, Ranges.get(i, j));
            }
            
            
            // Reads data for HM model
            
            u = getClass().getResource("/data/bcfcaesar_HM_centres.csv");
            DataCentresReader = new RBFNNDataReader(u.openStream(),0,0,0);
            u = getClass().getResource("/data/bcfcaesar_HM_weights.csv");
            DataWeightsReader = new RBFNNDataReader(u.openStream(),0,0,0);
            u = getClass().getResource("/data/bcfcaesar_HM_subset.csv");
            DataSubsetReader = new RBFNNDataReader(u.openStream(),0,0,0);

            HMCentres = new Matrix(DataCentresReader.Data);
            HMWeights = new Matrix(DataWeightsReader.Data);
            DataSubset = DataSubsetReader.DataAsInt();

            BufMatrix = new Matrix(HMCentres.getRowDimension(), DataSubsetReader.dim_j);
            for (int i=0; i<HMCentres.getRowDimension(); i++)
                for (int j=0; j<DataSubsetReader.dim_j; j++)
                    BufMatrix.set(i, j, HMCentres.get(i, DataSubset[0][j]));
            HMCentres = BufMatrix;

            HMRanges = new Matrix(2,5);
            for (int i=0; i<2; i++)  {
                HMRanges.set(i, 0, Ranges.get(i, 0));
                HMRanges.set(i, 1, Ranges.get(i, 5));
                HMRanges.set(i, 2, Ranges.get(i, 3));
                HMRanges.set(i, 3, Ranges.get(i, 6));
                HMRanges.set(i, 4, Ranges.get(i, 7));
            }
               
        } catch (Exception ex) {
            throw new InitFailureException(
                "Unable to load model parameters from local file");
        }        
    }
    
    
    public boolean Calculate(double[] Descriptors) {

        PredictionHM = Descriptor.MISSING_VALUE;
        PredictionGA = Descriptor.MISSING_VALUE;
        Prediction = Descriptor.MISSING_VALUE;
        
        // Checks descriptors
        for (double d : Descriptors)
            if (d == Descriptor.MISSING_VALUE)
                return false;
        
        Matrix GADataSet = new Matrix(1,5);
        Matrix HMDataSet = new Matrix(1,5);
        
        // Descriptors original order:
        // 0: X0sol
        // 1: MATS5v
        // 2: GATS5v
        // 3: BEHp2
        // 4: AEige
        // 5: Cl-089
        // 6: MlogP
        // 7: SsCl
        
        // Builds GADataset
        GADataSet.set(0, 0, Descriptors[6]); // MLogP
        GADataSet.set(0, 1, Descriptors[0]); // X0sol
        GADataSet.set(0, 2, Descriptors[1]); // MATS5v
        GADataSet.set(0, 3, Descriptors[3]); // BEHp2
        GADataSet.set(0, 4, Descriptors[7]); // SsCl 

        // Builds HMDataset
        HMDataSet.set(0, 0, Descriptors[6]); // MLogP
        HMDataSet.set(0, 1, Descriptors[4]); // AEige
        HMDataSet.set(0, 2, Descriptors[3]); // BEHp2
        HMDataSet.set(0, 3, Descriptors[2]); // GATS5v
        HMDataSet.set(0, 4, Descriptors[5]); // Cl-089
        
        try {

            // Prediction for GA and HM models
            RBFNNModel Model = new RBFNNModel(HMDataSet, HMCentres, HMWeights, HMRanges, 1.5);
            PredictionHM = Model.Predict().get(0, 0);

            // Prediction with HM descriptors model
            Model = new RBFNNModel(GADataSet, GACentres, GAWeights, GARanges, 1.2);
            PredictionGA = Model.Predict().get(0,0);

            // Hybrid model

            double ModelMean = (getPredictionHM() + getPredictionGA()) / 2;  
            double ModelMin = getPredictionGA();
            if (getPredictionGA() >= getPredictionHM())
                ModelMin = getPredictionHM();

            if (ModelMean>2.410)
                Prediction = (1.052 * ModelMean) - 0.065;
            else if  ( (ModelMean>1.355) && (ModelMean<=2.410) )
                Prediction = (0.996 * ModelMin) + 0.042;
            else 
                Prediction = (0.936 * ModelMean) - 0.123;

        } catch (Throwable e) {
            return false;
        }

        return true;
    }

    /**
     * @return the Prediction
     */
    public double getPrediction() {
        return Prediction;
    }

    /**
     * @return the PredictionHM
     */
    public double getPredictionHM() {
        return PredictionHM;
    }

    /**
     * @return the PredictionGA
     */
    public double getPredictionGA() {
        return PredictionGA;
    }
    
}
