package insilico.persistence_quantative_water_irfmn.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for ionization potential. Values are taken from Dragon 7.0
 * atom table.<p>
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class IonizationPotential {


    /**
     * Calculate ionization potential for all atoms in the given molecule, and returns
     * them in a single array.
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetIonizationPotential(at.getSymbol());
        }

        return w;
    }

    /**
     * Calculate ionization potential for all atoms in the given molecule, and returns
     * them in a single array. Weights normalized on C value.
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeightsNormalized(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];
        double wC = GetIonizationPotential("C");

        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            w[i] = GetIonizationPotential(at.getSymbol()) / wC;
        }

        return w;
    }

    /**
     * Returns the ionization potential for the given atom type.
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static double GetIonizationPotential(String AtomType) {

        if (AtomType.compareTo("H")==0) {
            return 13.5984;
        }
        if (AtomType.compareTo("B")==0) {
            return 8.298;
        }
        if (AtomType.compareTo("C")==0) {
            return 11.2603;
        }
        if (AtomType.compareTo("N")==0) {
            return 14.5341;
        }
        if (AtomType.compareTo("O")==0) {
            return 13.6181;
        }
        if (AtomType.compareTo("F")==0) {
            return 17.4228;
        }
        if (AtomType.compareTo("Al")==0) {
            return 5.9858;
        }
        if (AtomType.compareTo("Si")==0) {
            return 8.1517;
        }
        if (AtomType.compareTo("P")==0) {
            return 10.4867;
        }
        if (AtomType.compareTo("S")==0) {
            return 10.36;
        }
        if (AtomType.compareTo("Cl")==0) {
            return 12.9676;
        }
        if (AtomType.compareTo("Fe")==0) {
            return 7.9024;
        }
        if (AtomType.compareTo("Co")==0) {
            return 7.881;
        }
        if (AtomType.compareTo("Ni")==0) {
            return 7.6398;
        }
        if (AtomType.compareTo("Cu")==0) {
            return 7.7264;
        }
        if (AtomType.compareTo("Zn")==0) {
            return 9.3942;
        }
        if (AtomType.compareTo("Br")==0) {
            return 11.8138;
        }
        if (AtomType.compareTo("Sn")==0) {
            return 7.3438;
        }
        if (AtomType.compareTo("I")==0) {
            return 10.4513;
        }
        if (AtomType.compareTo("Gd")==0) {
            return 6.1498;
        }
        if (AtomType.compareTo("Cr")==0) {
            return 6.7665;
        }
        if (AtomType.compareTo("Mn")==0) {
            return 7.434;
        }
        if (AtomType.compareTo("Mo")==0) {
            return 7.0924;
        }
        if (AtomType.compareTo("Ag")==0) {
            return 7.5762;
        }
        if (AtomType.compareTo("Cd")==0) {
            return 8.9938;
        }
        if (AtomType.compareTo("Pt")==0) {
            return 8.9588;
        }
        if (AtomType.compareTo("Au")==0) {
            return 9.2255;
        }
        if (AtomType.compareTo("Hg")==0) {
            return 10.4375;
        }
        if (AtomType.compareTo("Se")==0) {
            return 9.7524;
        }
        if (AtomType.compareTo("Te")==0) {
            return 9.0096;
        }
        if (AtomType.compareTo("As")==0) {
            return 9.8152;
        }
        if (AtomType.compareTo("Ga")==0) {
            return 5.9993;
        }
        if (AtomType.compareTo("In")==0) {
            return 5.7864;
        }
        if (AtomType.compareTo("Tl")==0) {
            return 6.1082;
        }
        if (AtomType.compareTo("Pb")==0) {
            return 7.4166;
        }
        if (AtomType.compareTo("Ge")==0) {
            return 7.9;
        }
        if (AtomType.compareTo("Sb")==0) {
            return 8.6084;
        }
        if (AtomType.compareTo("Bi")==0) {
            return 7.2855;
        }

        return Descriptor.MISSING_VALUE;
    }





}
