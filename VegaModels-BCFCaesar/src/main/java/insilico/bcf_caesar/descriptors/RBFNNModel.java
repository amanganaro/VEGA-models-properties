package insilico.bcf_caesar.descriptors;

import Jama.Matrix;
import insilico.core.exception.GenericFailureException;


/**
 * Class ModelRBFNN implements a Radial Basis Function (RBF) Neural
 * Networks, providing a prediction upon a given dataset after having
 * received the neurons (centres) matrix and their weights.
 *
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class RBFNNModel {

    private Matrix Centres;         // Matrix (n x m) of NN centres (neurons)
    private Matrix Weights;         // Matrix (m x 1) of neurons weights
    private Matrix Data;            // Matrix (p x n) of p objects for n variables
    private Matrix Ranges;          // Matrix (2 x n) of ranges used to normalize
                                    // the n variables of Data
    private double ScaleFactor;     // Factor for scaling distances in RBF


    /**
     * Constructor of the class ModelRBFNN
     * 
     * @param Dataset   Matrix (p x n) of p objects for n variables
     * @param Centres   Matrix (n x m) of NN centres (neurons)
     * @param Weights   Matrix (m x 1) of neurons weights
     * @param NormalizeRanges   Matrix (2 x n) of ranges used to normalize
     *                          the n variables of Data
     * @param ScaleFactor   Factor for scaling distances in RBF
     */
    public RBFNNModel(Matrix Dataset, Matrix Centres, Matrix Weights,
                      Matrix NormalizeRanges, double ScaleFactor)
            throws GenericFailureException {

	// Sets object variables
        this.Centres = Centres;
        this.Weights = Weights;
        this.Data = Dataset;
        this.Ranges = NormalizeRanges;
        this.ScaleFactor = ScaleFactor;

        // Checks input matrices dimension
        if (Centres.getRowDimension()!=Data.getColumnDimension())
            throw new GenericFailureException(
                "Error in model: Neurons and dataset dimension not matching");
        if (Weights.getRowDimension()!=Centres.getColumnDimension())
            throw new GenericFailureException(
                "Error in model: Neurons and weights dimension not matching");
        if (NormalizeRanges.getColumnDimension()!=Data.getColumnDimension())
            throw new GenericFailureException(
                "Error in model: Dataset and ranges dimension not matching");
        if (NormalizeRanges.getRowDimension()!=2)
            throw new GenericFailureException(
                "Error in model: Wrong ranges dimension");

        // Normalizes data
        NormalizeData();

        // Transposes original dataset
        Data = Data.transpose();
	
    }


    /**
     * Normalizes the Data matrix using the ranges contained in Ranges
     */
    private void NormalizeData() {
        for (int j=0; j<Data.getColumnDimension(); j++)
            for (int i=0; i<Data.getRowDimension(); i++)
                Data.set(i, j, (-1 + 2 * (Data.get(i,j) - Ranges.get(1,j))/
		  (Ranges.get(0,j) - Ranges.get(1,j)) ));
    }


    /**
     * Builds the RBF design matrix on the Data matrix
     * 
     * @return  design matrix H (p x m)
     */
    private Matrix GetRBFDesign() {
		
        Matrix H;
        Matrix D;
        Matrix diag; //Vector diag;
        Matrix s; //Vector s;
        Matrix Buf;
        Matrix h;
        int n, p, m;
        int i, j, k, w;

        n = Data.getRowDimension();
        p = Data.getColumnDimension();
        m = Centres.getColumnDimension();

        // Create the H (p x m) result design matrix, filled with zeros 
        H = new Matrix(p, m);
        for (i=0; i<p; i++)
                for (j=0; j<m; j++)
                        H.set(i, j, 0);

        // Cycle for each centre (neuron)
        for (j=0; j<m; j++) {

                // Gets p difference vectors of each object from the
                // the j-th centre (neuron)
                D = new Matrix(n,p);
                for (k=0; k<p; k++)
                        for (w=0; w<n; w++)
                                //D.matrix[w][k] = Data.matrix[w][k] - Centres.matrix[w][j]; 
                                D.set(w, k, Data.get(w, k) - Centres.get(w, j)); 

                // Does metric calculation, obtaining the distance of each object
                // from the j-th centre (neuron), scaled with ScaleFactor
                // s = diag(D'*D) / scale^2
                diag = new Matrix(p, 1); //diag = new Vector(p);
                s = new Matrix(p,1);  //s = new Vector(p);
                Buf = new Matrix(p,p);
                Buf = (D.transpose()).times(D);
                for (int l=0; l<p; l++)
                    diag.set(l, 0, Buf.get(l, l));
                s = diag.times( 1 / Math.pow(ScaleFactor,2) );
//                diag = ( (D.transpose()).mul(D) ).getVectorFromDiagonal();
//                s = diag.mul( 1 / Math.pow(ScaleFactor,2) );

                // Applies radial basis function to distances
                h = new Matrix(1,p);
                for (k=0; k<p; k++)
                        h.set(0, k, Math.exp(s.get(k, 0) * (-1)));
//                h = new Matrix(1,s.size);
//                for (k=0; k<s.size; k++)
//                        h.matrix[0][k] = Math.exp(s.vector[k] * (-1));

                // Puts resulting distances vector h into the final H array
                h = h.transpose();
                for (k=0; k<p; k++)
//                    H.matrix[k][j] = h.matrix[k][0]; 
                    H.set(k, j, h.get(k, 0)); 
        }

        // Returns RBF design matrix
        return H;
    }


    /**
     * Makes prediction on Data using the NN built with the given centres
     *
     * @return  results matrix (p x 1)
     */
    public Matrix Predict() {

        Matrix H;
	Matrix Predictions;

	// Gets the design matrix from the input data, centre positions
	// and radius with Gaussian RBF
	H = GetRBFDesign();

	// Calculate predictions as (H * Weights)
//	Predictions = H.mul(Weights);
	Predictions = H.times(Weights);

	return Predictions;
    }

}
