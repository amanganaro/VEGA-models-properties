package insilico.fish_nic_tk.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for quantum number.
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class QuantumNumber {


    /**
     * Calculate quantum number for all atoms in the given molecule,
     * and returns them in a single array.
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static int[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        int[] w = new int[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetQuantumNumber(at.getSymbol());
        }

        return w;
    }

    /**
     * Returns the quantum number for the given atom type.
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static int GetQuantumNumber(String AtomType) {

        if (AtomType.compareTo("H")==0) {
            return 1;
        }
        if (AtomType.compareTo("B")==0) {
            return 5;
        }
        if (AtomType.compareTo("C")==0) {
            return 2;
        }
        if (AtomType.compareTo("N")==0) {
            return 2;
        }
        if (AtomType.compareTo("O")==0) {
            return 2;
        }
        if (AtomType.compareTo("F")==0) {
            return 2;
        }
        if (AtomType.compareTo("Al")==0) {
            return 3;
        }
        if (AtomType.compareTo("Si")==0) {
            return 3;
        }
        if (AtomType.compareTo("P")==0) {
            return 3;
        }
        if (AtomType.compareTo("S")==0) {
            return 3;
        }
        if (AtomType.compareTo("Cl")==0) {
            return 3;
        }
        if (AtomType.compareTo("Fe")==0) {
            return 4;
        }
        if (AtomType.compareTo("Co")==0) {
            return 4;
        }
        if (AtomType.compareTo("Ni")==0) {
            return 4;
        }
        if (AtomType.compareTo("Cu")==0) {
            return 4;
        }
        if (AtomType.compareTo("Zn")==0) {
            return 4;
        }
        if (AtomType.compareTo("Br")==0) {
            return 4;
        }
        if (AtomType.compareTo("Cr")==0) {
            return 4;
        }
        if (AtomType.compareTo("Mn")==0) {
            return 4;
        }
        if (AtomType.compareTo("Mo")==0) {
            return 5;
        }
        if (AtomType.compareTo("Se")==0) {
            return 4;
        }
        if (AtomType.compareTo("As")==0) {
            return 4;
        }
        if (AtomType.compareTo("Ga")==0) {
            return 4;
        }
        if (AtomType.compareTo("Ge")==0) {
            return 4;
        }
        if (AtomType.compareTo("Sn")==0) {
            return 5;
        }
        if (AtomType.compareTo("I")==0) {
            return 5;
        }

        return Descriptor.MISSING_VALUE;
    }

}
