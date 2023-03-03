package insilico.eye_irritation.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.weights.iWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsQuantumNumber;
import insilico.core.localization.StringSelectorCore;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *  Intrinsic State.
 *  Calculated for each atom as I = ( (2/L)^2 * delta_v + 1 ) / delta
 *  where:
 *  L is the principal quantum number of the atom
 *  delta_v is the valence vertex degree of the atom
 *  delta is the number of sigma electrons (simple vertex degree) of the atom
 *
 */
public class WeightsIState implements iWeight {

    private static final String SYMBOL = "s";
    private static final String NAME = StringSelectorCore.getString("descriptors_istate_name");

    private final static Object[][] EL_NUMBER = {
            {"H",1},
            {"He",1},
            {"Li",1},
            {"Be",2},
            {"B",3},
            {"C",4},
            {"N",5},
            {"O",6},
            {"F",7},
            {"Ne",8},
            {"Na",1},
            {"Mg",2},
            {"Al",3},
            {"Si",4},
            {"P",5},
            {"S",6},
            {"Cl",7},
            {"Ar",8},
            {"K",1},
            {"Ca",2},
            {"Sc",3},
            {"Ti",4},
            {"V",5},
            {"Cr",6},
            {"Mn",7},
            {"Fe",8},
            {"Co",9},
            {"Ni",10},
            {"Cu",11},
            {"Zn",12},
            {"Ga",3},
            {"Ge",4},
            {"As",5},
            {"Se",6},
            {"Br",7},
            {"Kr",8},
            {"Rb",1},
            {"Sr",2},
            {"Y",3},
            {"Zr",4},
            {"Nb",5},
            {"Mo",6},
            {"Tc",7},
            {"Ru",8},
            {"Rh",9},
            {"Pd",10},
            {"Ag",11},
            {"Cd",12},
            {"In",3},
            {"Sn",4},
            {"Sb",5},
            {"Te",6},
            {"I",7},
            {"Xe",8},
            {"Cs",1},
            {"Ba",2},
            {"La",3},
            {"Ce",4},
            {"Pr",5},
            {"Nd",6},
            {"Pm",7},
            {"Sm",8},
            {"Eu",9},
            {"Gd",10},
            {"Tb",11},
            {"Dy",12},
            {"Ho",13},
            {"Er",14},
            {"Tm",15},
            {"Yb",16},
            {"Lu",17},
            {"Hf",4},
            {"Ta",5},
            {"W",6},
            {"Re",7},
            {"Os",8},
            {"Ir",9},
            {"Pt",10},
            {"Au",11},
            {"Hg",12},
            {"Tl",3},
            {"Pb",4},
            {"Bi",5},
            {"Po",6},
            {"At",7},
            {"Rn",8},
            {"Fr",1},
            {"Ra",2},
            {"Ac",3},
            {"Th",4},
            {"Pa",5},
            {"U",6},
            {"Np",7},
            {"Pu",8},
            {"Am",9},
            {"Cm",10},
            {"Bk",11},
            {"Cf",12},
            {"Es",13},
            {"Fm",14},
            {"Md",15},
            {"No",16},
            {"Lr",17},
    };


    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public String getSymbol() {
        return SYMBOL;
    }


    public double[] getWeights(IAtomContainer mol, boolean HasExplicitH)  {

        int nSK = mol.getAtomCount();
        double[] w = new double[nSK];

        for (int i=0; i<nSK; i++)
            w[i] = getWeight(mol, i, HasExplicitH);

        return w;
    }


    public double getWeight(IAtomContainer mol, int atomIndex, boolean HasExplicitH) {

        IAtom a = mol.getAtom(atomIndex);

        //// Quantum number / period of the atom
        int L = (new WeightsQuantumNumber()).getWeight(a.getSymbol()) ;


        //// Electron number
        int nElectrons = Descriptor.MISSING_VALUE;
        for (Object[] buf : EL_NUMBER) {
            String symbol = (String) buf[0];
            if (symbol.equalsIgnoreCase(a.getSymbol())) {
                nElectrons = (int) buf[1];
                break;
            }
        }
        if (nElectrons == Descriptor.MISSING_VALUE)
            return Descriptor.MISSING_VALUE;


        //// Vertex degree of the atom
        double D = 0;
        if (HasExplicitH) {
            int nH = 0;
            for (IAtom conn : mol.getConnectedAtomsList(a))
                if (conn.getAtomicNumber() == 1)
                    nH++;
            D = mol.getConnectedBondsCount(a) - nH;
        } else {
            D = mol.getConnectedBondsCount(a);
        }


        //// Valence vertex degree of the atom
        int nH = 0;
        if (HasExplicitH) {
            for (IAtom conn : mol.getConnectedAtomsList(a))
                if (conn.getAtomicNumber() == 1)
                    nH++;
        } else {
            nH = a.getImplicitHydrogenCount();
        }
        double DeltaV = nElectrons - nH;



        //// I-State
        double IS = (Math.pow(2.0 / L, 2.0) * DeltaV + 1) / D;
        if ( (Double.isInfinite(IS)) || (Double.isNaN(IS)) )
            IS = Descriptor.MISSING_VALUE;
        return IS;
    }
}
