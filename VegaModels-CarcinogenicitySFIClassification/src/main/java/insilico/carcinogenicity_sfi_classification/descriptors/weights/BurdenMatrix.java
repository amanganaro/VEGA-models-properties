package insilico.carcinogenicity_sfi_classification.descriptors.weights;

import insilico.core.molecule.matrix.ConnectionAugMatrix;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Calculates Burden matrix.<p>
 * Each element on the diagonal has the value of Z for the i-th atom
 * Each element (i,j) outside the diagonal has the value of the square root
 * of the bond order between i and j, or 0.001 if they are not adjacent.
 * Aromatic bond is coded as 1.5
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class BurdenMatrix {
    /**
     * Calculates Burden matrix.<p>
     * Each element on the diagonal has the value of Z for the i-th atom
     * Each element (i,j) outside the diagonal has the value of the square root
     * of the bond order between i and j, or 0.001 if they are not adjacent.
     * Aromatic bond is coded as 1.5
     *
     * @param molecule source CDK Molecule
     * @return augmented connection matrix
     */
    static public double[][] getMatrix(IAtomContainer molecule) {

        int nSK = molecule.getAtomCount();
        double[][] matrix = new double[nSK][nSK];

        double[][] ConnMat = ConnectionAugMatrix.getMatrix(molecule);

        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {

                if (i==j) {
                    matrix[i][j] = ConnMat[i][j];
                } else {
                    if (ConnMat[i][j] > 0) {
                        matrix[i][j] = Math.sqrt(ConnMat[i][j]);
                    } else {
                        matrix[i][j] = 0.001;
                    }
                }

            }

        return matrix;
    }
}
