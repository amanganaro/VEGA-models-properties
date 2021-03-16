package insilico.daphnia_epa.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.molecule.tools.AtomicNumber;
import insilico.core.molecule.tools.Manipulator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Weighting scheme for valence vertex degree.
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ValenceVertexDegree {

    /**
     * Calculate the valence vertex degree for all atoms in the given molecule,
     * and returns them in a single array.
     *
     * @param mol molecule to be processed
     * @return array of weights for all molecule's atoms
     */
    public static double[] getWeights(IAtomContainer mol) {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];
        AtomicNumber ZFinder;
        try {
            ZFinder = new AtomicNumber();
        } catch (Exception e) {
            for (int i=0; i<nSK; i++)
                w[i] = -999;
            return w;
        }

        for (int i=0; i<nSK; i++) {
            IAtom at =  mol.getAtom(i);
            int Z = ZFinder.GetAtomicNumber(at.getSymbol());
            int Zv = GetValenceElectronsNumber(at.getSymbol());
            int h = Manipulator.CountImplicitHydrogens(at);
            int ch = at.getFormalCharge();
            if (Zv == -999)
                w[i] = -999;
            else
                w[i] = (double)(Zv - h - ch) / (double)(Z - Zv - 1.00);
        }

        return w;
    }

    /**
     * Returns the valence electron number for the given atom type.
     *
     * @param AtomType atomtype symbol to be processed
     * @return the value of the weight
     */
    public static int GetValenceElectronsNumber(String AtomType) {

        if (AtomType.compareToIgnoreCase("H")==0) {
            return 1;
        }
        if (AtomType.compareToIgnoreCase("B")==0) {
            return 3;
        }
        if (AtomType.compareToIgnoreCase("C")==0) {
            return 4;
        }
        if (AtomType.compareToIgnoreCase("N")==0) {
            return 5;
        }
        if (AtomType.compareToIgnoreCase("O")==0) {
            return 6;
        }
        if (AtomType.compareToIgnoreCase("F")==0) {
            return 7;
        }
        if (AtomType.compareToIgnoreCase("Al")==0) {
            return 3;
        }
        if (AtomType.compareToIgnoreCase("Si")==0) {
            return 4;
        }
        if (AtomType.compareToIgnoreCase("P")==0) {
            return 5;
        }
        if (AtomType.compareToIgnoreCase("S")==0) {
            return 6;
        }
        if (AtomType.compareToIgnoreCase("Cl")==0) {
            return 7;
        }
        if (AtomType.compareToIgnoreCase("Fe")==0) {
            return 8;
        }
        if (AtomType.compareToIgnoreCase("Co")==0) {
            return 9;
        }
        if (AtomType.compareToIgnoreCase("Ni")==0) {
            return 10;
        }
        if (AtomType.compareToIgnoreCase("Cu")==0) {
            return 11;
        }
        if (AtomType.compareToIgnoreCase("Zn")==0) {
            return 12;
        }
        if (AtomType.compareToIgnoreCase("Br")==0) {
            return 7;
        }
        if (AtomType.compareToIgnoreCase("Cr")==0) {
            return 6;
        }
        if (AtomType.compareToIgnoreCase("Mn")==0) {
            return 7;
        }
        if (AtomType.compareToIgnoreCase("Mo")==0) {
            return 6;
        }
        if (AtomType.compareToIgnoreCase("Se")==0) {
            return 6;
        }
        if (AtomType.compareToIgnoreCase("As")==0) {
            return 5;
        }
        if (AtomType.compareToIgnoreCase("Ga")==0) {
            return 3;
        }
        if (AtomType.compareToIgnoreCase("Ge")==0) {
            return 4;
        }
        if (AtomType.compareToIgnoreCase("Sn")==0) {
            return 4;
        }
        if (AtomType.compareToIgnoreCase("I")==0) {
            return 7;
        }

        return Descriptor.MISSING_VALUE;
    }
}
