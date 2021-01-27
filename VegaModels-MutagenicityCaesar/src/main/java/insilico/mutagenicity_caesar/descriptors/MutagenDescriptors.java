package insilico.mutagenicity_caesar.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * Class that implements the 21 fragment descriptors needed by the Mutagen model:
 * 0 - SsCH3
 * 1 - SdCH2
 * 2 - SssCH2
 * 3 - SdsCH
 * 4 - SaaCH
 * 5 - SsssCH
 * 6 - SdssC
 * 7 - SaasC
 * 8 - SaaaC
 * 9 - SssssC
 * 10 - SsNH2
 * 11 - StN
 * 12 - SdsN
 * 13 - SaaN
 * 14 - SsssN
 * 15 - SsaaN
 * 16 - SsOH
 * 17 - SdO
 * 18 - SssO
 * 19 - SaaO
 * 20 - SHCHnX
 * 
 * @author amanganaro
 */
public class MutagenDescriptors {

    /**
     * Calculates the descriptors for the given Molecule CurMol, returns
     * them into a double vector.
     * 
     * @param mol
     * @return
     */
    public static double[] Calculate(InsilicoMolecule mol) {
    
        double[] descriptors = new double[21];
        for (int i=0; i<21; i++) descriptors[i]=0;
        
        try {
            
            IAtomContainer Mol = mol.GetStructure();

            //
            // Calculates descriptors
            //

            int nAtoms = Mol.getAtomCount();
            int[] nSingle=new int[nAtoms], nDouble=new int[nAtoms], 
                nTriple=new int[nAtoms], nArom=new int[nAtoms], nH=new int[nAtoms],
                atType=new int[nAtoms], atArom=new int[nAtoms],
                CarbonHaloGroup=new int[nAtoms];

            for (int i=0; i<nAtoms; i++) {
                nSingle[i] = 0;
                nDouble[i] = 0;
                nTriple[i] = 0;
                nH[i] = 0;
                nArom[i] = 0;
                atArom[i] = 0;
                atType[i] = SymbolToAtNum(Mol.getAtom(i).getSymbol());
                CarbonHaloGroup[i] = 0;
            }

            for (int i=0; i<Mol.getBondCount(); i++) {

                int bondOrder = 0;
                IBond CurBond =  Mol.getBond(i);

                int atomA = SymbolToAtNum((CurBond.getAtom(0)).getSymbol());
                int atomB = SymbolToAtNum((CurBond.getAtom(1)).getSymbol());
                int numA = Mol.indexOf(CurBond.getAtom(0));
                int numB = Mol.indexOf(CurBond.getAtom(1));

                if (CurBond.getFlag(CDKConstants.ISAROMATIC))
                    bondOrder = 9;
                else {
                    if (CurBond.getOrder()== Order.SINGLE) 
                        bondOrder = 1;
                    else if (CurBond.getOrder()== Order.DOUBLE) 
                        bondOrder = 2;
                    else if (CurBond.getOrder()== Order.TRIPLE) 
                        bondOrder = 3;
                    else if (CurBond.getOrder()== Order.QUADRUPLE) 
                        bondOrder = 4;

                    if ( ( (CurBond.getAtom(0).getFormalCharge()==1)&&
                       (CurBond.getAtom(1).getFormalCharge()==-1) ) || 
                       ( (CurBond.getAtom(1).getFormalCharge()==1)&&
                       (CurBond.getAtom(0).getFormalCharge()==-1) )){
                        bondOrder++;
                        CurBond.getAtom(0).setFormalCharge(0);
                        CurBond.getAtom(1).setFormalCharge(0);
                    }
                }

                if ((atomA>=6)&&(atomA<=8)) {
                    nH[numA]=CurBond.getAtom(0).getImplicitHydrogenCount();
                    switch (bondOrder) {
                        case 1: nSingle[numA]++; break;
                        case 2: nDouble[numA]++; break;
                        case 3: nTriple[numA]++; break;
                        case 9: nArom[numA]++; break;
                    }
                }

                if ((atomB>=6)&&(atomB<=8)) {
                    nH[numB]=CurBond.getAtom(1).getImplicitHydrogenCount();
                    switch (bondOrder) {
                        case 1: nSingle[numB]++; break;
                        case 2: nDouble[numB]++; break;
                        case 3: nTriple[numB]++; break;
                        case 9: nArom[numB]++; break;
                    }
                }

                if ((atomA==6)&&(atomB==9)) 
                    CarbonHaloGroup[numA]++;

                if ((atomB==6)&&(atomA==9)) 
                    CarbonHaloGroup[numB]++;

            }


            for (int i=0; i<nAtoms; i++) {

                switch (atType[i]) {
                    
                    case 6: // C atoms

                        // 0 - SsCH3
                        if ((nH[i]==3) && (nSingle[i]==1))
                            descriptors[0]++;

                        // 1 - SdCH2
                        if ((nH[i]==2) && (nDouble[i]==1))
                            descriptors[1]++;

                        // 2 - SssCH2
                        if ((nH[i]==2) && (nSingle[i]==2))
                            descriptors[2]++;

                        // 3 - SdsCH
                        if ((nH[i]==1) && (nSingle[i]==1) && (nDouble[i]==1))
                            descriptors[3]++;

                        // 4 - SaaCH
                        if ((nH[i]==1) && (nArom[i]==2))
                            descriptors[4]++;

                        // 5 - SsssCH
                        if ((nH[i]==1) && (nSingle[i]==3))
                            descriptors[5]++;

                        // 6 - SdssC
                        if ((nH[i]==0) && (nSingle[i]==2) && (nDouble[i]==1))
                            descriptors[6]++;

                        // 7 - SaasC
                        if ((nH[i]==0) && (nSingle[i]==1) && (nArom[i]==2))
                            descriptors[7]++;

                        // 8 - SaaaC
                        if ((nH[i]==0) && (nArom[i]==3))
                            descriptors[8]++;

                        // 9 - SssssC
                        if ((nH[i]==0) && (nSingle[i]==4))
                            descriptors[9]++;

                        // 20 - SHCHnX
                        if ( ((nH[i]==1)||(nH[i]==2)) && (CarbonHaloGroup[i]>0))
                            descriptors[20]++;

                        break;
                        
                    case 7: // N atoms

                        // 10 - SsNH2
                        if ((nH[i]==2) && (nSingle[i]==1) && (nDouble[i]==0) && (nArom[i]==0))
                            descriptors[10]++;

                        // 11 - StN
                        if ((nH[i]==0) && (nSingle[i]==0) && (nDouble[i]==0) && (nArom[i]==0) &&
                                (nTriple[i]==1))
                            descriptors[11]++;

                        // 12 - SdsN
                        if ((nH[i]==0) && (nSingle[i]==1) && (nDouble[i]==1) && (nArom[i]==0))
                            descriptors[12]++;

                        // 13 - SaaN
                        if ((nH[i]==0) && (nSingle[i]==0) && (nDouble[i]==0) && (nArom[i]==2))
                            descriptors[13]++;

                        // 14 - SsssN
                        if ((nH[i]==0) && (nSingle[i]==3) && (nDouble[i]==0) && (nArom[i]==0))
                            descriptors[14]++;

                        // 15 - SsaaN
                        if ((nH[i]==0) && (nSingle[i]==1) && (nDouble[i]==0) && (nArom[i]==2))
                            descriptors[15]++;

                        break;
                        
                    case 8: // O atoms

                        // 16 - SsOH
                        if ((nH[i]==1) && (nSingle[i]==1))
                            descriptors[16]++;

                        // 17 - SdO
                        if ((nH[i]==0) && (nDouble[i]==1))
                            descriptors[17]++;

                        // 18 - SssO
                        if ((nH[i]==0) && (nSingle[i]==2))
                            descriptors[18]++;

                        // 19 - SaaO
                        if ((nH[i]==0) && (nArom[i]==2))
                            descriptors[19]++;

                        break;
                }

            }

        } catch (Exception e) {
            for (int i=0; i<24; i++)
                descriptors[i] = Descriptor.MISSING_VALUE;
        }
         
        
        return descriptors;
        
    }


    
    private static int SymbolToAtNum(String Symbol) {
        
        int ret = 0;
        
        if (Symbol.compareToIgnoreCase("C")==0)
            ret = 6;
        if (Symbol.compareToIgnoreCase("N")==0)
            ret = 7;
        if (Symbol.compareToIgnoreCase("O")==0)
            ret = 8;
        if (Symbol.compareToIgnoreCase("F")==0)
            ret = 9; // Halogen
        if (Symbol.compareToIgnoreCase("Cl")==0)
            ret = 9; // Halogen
        
        return ret;
    }
    
}
