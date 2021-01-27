package insilico.fathead_epa.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for electronegativity. Values are taken from Dragon 5.5
 * documentation. NOTE: values are scaled on the reference value (Carbon)
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class Electronegativity {

    /**
     * Calculate electronegativity weights for all atoms in the given molecule,
     * and returns them in a single array. NOTE: values are scaled on the reference
     * value (Carbon).
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetElectronegativity(at.getSymbol());
        }

        return w;
    }


    /**
     * Returns the electronegativity weight for the given atom type. NOTE: values are
     * scaled on the reference value (Carbon).
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static double GetElectronegativity(String AtomType) {

        if (AtomType.compareTo("H")==0) {
            return 0.94;
        }
        if (AtomType.compareTo("B")==0) {
            return 0.83;
        }
        if (AtomType.compareTo("C")==0) {
            return 1;
        }
        if (AtomType.compareTo("N")==0) {
            return 1.16;
        }
        if (AtomType.compareTo("O")==0) {
            return 1.33;
        }
        if (AtomType.compareTo("F")==0) {
            return 1.451;
        }
        if (AtomType.compareTo("Al")==0) {
            return 0.62;
        }
        if (AtomType.compareTo("Si")==0) {
            return 0.78;
        }
        if (AtomType.compareTo("P")==0) {
            return 0.92;
        }
        if (AtomType.compareTo("S")==0) {
            return 1.08;
        }
        if (AtomType.compareTo("Cl")==0) {
            return 1.27;
        }
        if (AtomType.compareTo("Fe")==0) {
            return 0.80;
        }
        if (AtomType.compareTo("Co")==0) {
            return 0.93;
        }
        if (AtomType.compareTo("Ni")==0) {
            return 0.71;
        }
        if (AtomType.compareTo("Cu")==0) {
            return 0.72;
        }
        if (AtomType.compareTo("Zn")==0) {
            return 0.81;
        }
        if (AtomType.compareTo("Br")==0) {
            return 1.17;
        }
        if (AtomType.compareTo("Sn")==0) {
            return 0.84;
        }
        if (AtomType.compareTo("I")==0) {
            return 1.01;
        }
        if (AtomType.compareTo("Gd")==0) {
            return 0.73;
        }
        if (AtomType.compareTo("Cr")==0) {
            return 0.60;
        }
        if (AtomType.compareTo("Mn")==0) {
            return 0.80;
        }
        if (AtomType.compareTo("Mo")==0) {
            return 0.42;
        }
        if (AtomType.compareTo("Ag")==0) {
            return 0.67;
        }
        if (AtomType.compareTo("Cd")==0) {
            return 0.72;
        }
        if (AtomType.compareTo("Pt")==0) {
            return 0.83;
        }
        if (AtomType.compareTo("Au")==0) {
            return 0.92;
        }
        if (AtomType.compareTo("Hg")==0) {
            return 0.80;
        }
        if (AtomType.compareTo("Se")==0) {
            return 1.09;
        }
        if (AtomType.compareTo("Te")==0) {
            return 0.95;
        }
        if (AtomType.compareTo("As")==0) {
            return 1.03;
        }
        if (AtomType.compareTo("Ga")==0) {
            return 0.88;
        }
        if (AtomType.compareTo("In")==0) {
            return 0.78;
        }
        if (AtomType.compareTo("Tl")==0) {
            return 0.82;
        }
        if (AtomType.compareTo("Pb")==0) {
            return 0.83;
        }
        if (AtomType.compareTo("Ge")==0) {
            return 0.95;
        }
        if (AtomType.compareTo("Sb")==0) {
            return 0.89;
        }
        if (AtomType.compareTo("Bi")==0) {
            return 0.85;
        }

        return Descriptor.MISSING_VALUE;
    }
}
