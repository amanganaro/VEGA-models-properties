package insilico.pxr_up.descriptors.weights;

import insilico.core.molecule.tools.Manipulator;
import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for vertex degree. 
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class VertexDegree {

    
    /**
     * Calculate vertex degree for all atoms in the given molecule, 
     * and returns them in a single array.
     * 
     * @param mol molecule to be processed
     * @param HDepleted true if calculation should be on H depleted molecule,
     * otherwise also (implicit) hydrogen atoms are counted in the degree
     * @return array of weights for all molecule's atoms
     */
    public static int[] getWeights(IAtomContainer mol, boolean HDepleted) {
        
        int nSK = mol.getAtomCount();
        int[] w = new int[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at =  mol.getAtom(i);
            w[i] = mol.getConnectedAtomsCount(at);
            if (!(HDepleted))
                w[i] += Manipulator.CountImplicitHydrogens(at);
        }
        
        return w;
    }
}
