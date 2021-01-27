package insilico.daphnia_demetra.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for polarizability. Values are taken from Dragon 5.5
 * documentation. NOTE: values are scaled on the reference value (Carbon)
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class Polarizability {


    /**
     * Calculate polarizability weights for all atoms in the given molecule,
     * and returns them in a single array. NOTE: values are scaled on the
     * reference value (Carbon).
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetPolarizability(at.getSymbol());
        }

        return w;
    }


    /**
     * Returns the polarizability weight for the given atom type. NOTE: values are
     * scaled on the reference value (Carbon).
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static double GetPolarizability(String AtomType) {

        if (AtomType.compareTo("H")==0) {
            return 0.38;
        }
        if (AtomType.compareTo("B")==0) {
            return 1.72;
        }
        if (AtomType.compareTo("C")==0) {
            return 1;
        }
        if (AtomType.compareTo("N")==0) {
            return 0.62;
        }
        if (AtomType.compareTo("O")==0) {
            return 0.45;
        }
        if (AtomType.compareTo("F")==0) {
            return 0.31;
        }
        if (AtomType.compareTo("Al")==0) {
            return 3.86;
        }
        if (AtomType.compareTo("Si")==0) {
            return 3.06;
        }
        if (AtomType.compareTo("P")==0) {
            return 2.06;
        }
        if (AtomType.compareTo("S")==0) {
            return 1.65;
        }
        if (AtomType.compareTo("Cl")==0) {
            return 1.24;
        }
        if (AtomType.compareTo("Fe")==0) {
            return 4.77;
        }
        if (AtomType.compareTo("Co")==0) {
            return 4.26;
        }
        if (AtomType.compareTo("Ni")==0) {
            return 3.86;
        }
        if (AtomType.compareTo("Cu")==0) {
            return 3.47;
        }
        if (AtomType.compareTo("Zn")==0) {
            return 4.03;
        }
        if (AtomType.compareTo("Br")==0) {
            return 1.73;
        }
        if (AtomType.compareTo("Sn")==0) {
            return 4.38;
        }
        if (AtomType.compareTo("I")==0) {
            return 3.04;
        }
        if (AtomType.compareTo("Gd")==0) {
            return 13.35;
        }
        if (AtomType.compareTo("Cr")==0) {
            return 6.59;
        }
        if (AtomType.compareTo("Mn")==0) {
            return 5.34;
        }
        if (AtomType.compareTo("Mo")==0) {
            return 7.27;
        }
        if (AtomType.compareTo("Ag")==0) {
            return 4.09;
        }
        if (AtomType.compareTo("Cd")==0) {
            return 4.09;
        }
        if (AtomType.compareTo("Pt")==0) {
            return 3.69;
        }
        if (AtomType.compareTo("Au")==0) {
            return 3.30;
        }
        if (AtomType.compareTo("Hg")==0) {
            return 3.24;
        }
        if (AtomType.compareTo("Se")==0) {
            return 2.14;
        }
        if (AtomType.compareTo("Te")==0) {
            return 3.13;
        }
        if (AtomType.compareTo("As")==0) {
            return 2.45;
        }
        if (AtomType.compareTo("Ga")==0) {
            return 4.61;
        }
        if (AtomType.compareTo("In")==0) {
            return 5.80;
        }
        if (AtomType.compareTo("Tl")==0) {
            return 4.32;
        }
        if (AtomType.compareTo("Pb")==0) {
            return 3.86;
        }
        if (AtomType.compareTo("Ge")==0) {
            return 3.45;
        }
        if (AtomType.compareTo("Sb")==0) {
            return 3.75;
        }
        if (AtomType.compareTo("Bi")==0) {
            return 4.2;
        }

        return Descriptor.MISSING_VALUE;
    }

}
