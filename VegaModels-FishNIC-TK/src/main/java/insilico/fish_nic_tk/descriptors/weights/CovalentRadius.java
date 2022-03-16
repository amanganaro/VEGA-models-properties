package insilico.fish_nic_tk.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for Covalent Radius. Values are taken from Dragon 6
 * documentation.
 *
 * @author Alberto Manganaro
 */
public class CovalentRadius {

    /**
     * Calculate Covalent Radius weights for all atoms in the given
     * molecule, and returns them in a single array.
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetCovRadius(at.getSymbol());
        }

        return w;
    }


    /**
     * Returns the Covalent Radius weight for the given atom type.
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static double GetCovRadius(String AtomType) {

        if (AtomType.compareTo("H")==0) {
            return 0.31;
        }
        if (AtomType.compareTo("He")==0) {
            return 0.28;
        }
        if (AtomType.compareTo("Li")==0) {
            return 1.28;
        }
        if (AtomType.compareTo("Be")==0) {
            return 0.96;
        }
        if (AtomType.compareTo("B")==0) {
            return 0.84;
        }
        if (AtomType.compareTo("C")==0) {
            return 0.76;
        }
        if (AtomType.compareTo("N")==0) {
            return 0.71;
        }
        if (AtomType.compareTo("O")==0) {
            return 0.66;
        }
        if (AtomType.compareTo("F")==0) {
            return 0.57;
        }
        if (AtomType.compareTo("Ne")==0) {
            return 0.58;
        }
        if (AtomType.compareTo("Na")==0) {
            return 1.66;
        }
        if (AtomType.compareTo("Mg")==0) {
            return 1.41;
        }
        if (AtomType.compareTo("Al")==0) {
            return 1.21;
        }
        if (AtomType.compareTo("Si")==0) {
            return 1.11;
        }
        if (AtomType.compareTo("P")==0) {
            return 1.07;
        }
        if (AtomType.compareTo("S")==0) {
            return 1.05;
        }
        if (AtomType.compareTo("Cl")==0) {
            return 1.02;
        }
        if (AtomType.compareTo("Ar")==0) {
            return 1.06;
        }
        if (AtomType.compareTo("K")==0) {
            return 2.03;
        }
        if (AtomType.compareTo("Ca")==0) {
            return 1.76;
        }
        if (AtomType.compareTo("Sc")==0) {
            return 1.7;
        }
        if (AtomType.compareTo("Ti")==0) {
            return 1.6;
        }
        if (AtomType.compareTo("V")==0) {
            return 1.53;
        }
        if (AtomType.compareTo("Cr")==0) {
            return 1.39;
        }
        if (AtomType.compareTo("Mn")==0) {
            return 1.39;
        }
        if (AtomType.compareTo("Fe")==0) {
            return 1.32;
        }
        if (AtomType.compareTo("Con")==0) return 1.26;
        if (AtomType.compareTo("Ni")==0) return 1.24;
        if (AtomType.compareTo("Cu")==0) return 1.32;
        if (AtomType.compareTo("Zn")==0) return 1.22;
        if (AtomType.compareTo("Ga")==0) return 1.22;
        if (AtomType.compareTo("Ge")==0) return 1.2;
        if (AtomType.compareTo("As")==0) return 1.19;
        if (AtomType.compareTo("Se")==0) return 1.2;
        if (AtomType.compareTo("Br")==0) return 1.2;
        if (AtomType.compareTo("Kr")==0) return 1.16;
        if (AtomType.compareTo("Rb")==0) return 2.2;
        if (AtomType.compareTo("Sr")==0) return 1.95;
        if (AtomType.compareTo("Y")==0) return 1.9;
        if (AtomType.compareTo("Zr")==0) return 1.75;
        if (AtomType.compareTo("Nb")==0) return 1.64;
        if (AtomType.compareTo("Mo")==0) return 1.54;
        if (AtomType.compareTo("Tc")==0) return 1.47;
        if (AtomType.compareTo("Ru")==0) return 1.46;
        if (AtomType.compareTo("Rh")==0) return 1.42;
        if (AtomType.compareTo("Pd")==0) return 1.39;
        if (AtomType.compareTo("Ag")==0) return 1.45;
        if (AtomType.compareTo("Cd")==0) return 1.44;
        if (AtomType.compareTo("In")==0) return 1.42;
        if (AtomType.compareTo("Sn")==0) return 1.39;
        if (AtomType.compareTo("Sb")==0) return 1.39;
        if (AtomType.compareTo("Te")==0) return 1.38;
        if (AtomType.compareTo("I")==0) return 1.39;
        if (AtomType.compareTo("Xe")==0) return 1.4;
        if (AtomType.compareTo("Cs")==0) return 2.44;
        if (AtomType.compareTo("Ba")==0) return 2.15;
        if (AtomType.compareTo("La")==0) return 2.07;
        if (AtomType.compareTo("Ce")==0) return 2.04;
        if (AtomType.compareTo("")==0) return 0.0;

        return Descriptor.MISSING_VALUE;
    }
}
