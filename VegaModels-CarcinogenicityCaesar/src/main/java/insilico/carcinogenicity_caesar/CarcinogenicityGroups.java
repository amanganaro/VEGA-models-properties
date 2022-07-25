package insilico.carcinogenicity_caesar;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.log4j.Log4j;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class CarcinogenicityGroups {

    // Functional groups
    private int nRNNOx;
    private int nPO4;
    
    // Atom centered fragments
    private int N_067;
    private int N_078;
    
    
    
    public CarcinogenicityGroups(InsilicoMolecule Mol) {
        
        // Gets matrices
        double[][] ConnMatrix;
        try {
            ConnMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            nRNNOx = Descriptor.MISSING_VALUE;
            nPO4 = Descriptor.MISSING_VALUE;
            N_067 = Descriptor.MISSING_VALUE;
            N_078 = Descriptor.MISSING_VALUE;
            return;
        }
        
        int nSK;
        try {
            nSK = Mol.GetStructure().getAtomCount();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            nRNNOx = Descriptor.MISSING_VALUE;
            nPO4 = Descriptor.MISSING_VALUE;
            N_067 = Descriptor.MISSING_VALUE;
            N_078 = Descriptor.MISSING_VALUE;
            return;
        }
        
        nRNNOx = 0;
        nPO4 = 0;
        N_067 = 0;
        N_078 = 0;
        
        for (int i=0; i<nSK; i++) {

            // Atom: C
            if (ConnMatrix[i][i] == 6) {
                
                
                continue;
            }
            
            
            // Atom: N
            if (ConnMatrix[i][i] == 7) {
                
                int VD=0, nH=0, Charge=0;
                int nSingle=0, nDouble=0, nTriple=0, nArom=0;
                int nElNeg=0, sElNeg=0, dElNeg=0, nO=0, dO=0, nAromLinked=0; 
                int sC=0, sNdblO=0, sOR=0, sRdblX=0;
                
                for (int j=0; j<nSK; j++) {
                    if (j==i) 
                        continue;
                    if (ConnMatrix[i][j]>0) {
                        VD++;
                        double Z = ConnMatrix[j][j];
                        double b = ConnMatrix[i][j];

                        if (b==1) nSingle++;
                        if (b==2) nDouble++;
                        if (b==3) nTriple++;
                        if (b==1.5) nArom++;
                        
                        // Checks for electronegative atoms
                        if ((Z==7)||(Z==8)||(Z==15)||(Z==16)||(Z==34)||
                            (Z==9)||(Z==17)||(Z==35)||(Z==53)) {
                            nElNeg++;
                            if (b==1) sElNeg++;
                            if (b==2) dElNeg++;
                            if (Z==8) {
                                nO++;
                                if (b==2) dO++;
                                
                                // Checks for RCO-N
                                if (b==1) {
                                    int OVD=0;
                                    for (int k=0; k<nSK; k++) {
                                        if (k==j) continue;
                                        if (ConnMatrix[j][k]>0) 
                                            OVD++;
                                    }
                                    if (OVD == 2)
                                        sOR++;
                                }
                            }
                        }
                        
                        // Checks carbons
                        if (Z==6) {
                            if (b==1) sC++;
                        }
                        
                        // Checks for N=O a group
                        if ((Z==7)&&(b==1)) {
                            int LocalVD = 1;
                            boolean dblO = false;
                            for (int k=0; k<nSK; k++) {
                                if ((k==j)||(k==i)) 
                                    continue;
                                if (ConnMatrix[j][k]>0) {
                                    LocalVD++;
                                    if (ConnMatrix[j][k]==2)
                                        if (ConnMatrix[k][k]==8)
                                            dblO = true;
                                }
                            }
                            if ((LocalVD == 2) && (dblO))
                                sNdblO++;
                        }
                        
                        // Checks for N-R=X
                        if (b==1) {
                            for (int k=0; k<nSK; k++) {
                                if ((k==j)||(k==i)) continue;
                                if (ConnMatrix[j][k] == 2) {
                                    double LocalZ = ConnMatrix[k][k];
                                    if ((LocalZ==7)||(LocalZ==8)||(LocalZ==15)||
                                        (LocalZ==16)||(LocalZ==34)||(LocalZ==9)||
                                        (LocalZ==17)||(LocalZ==35)||(LocalZ==53))
                                        sRdblX++;
                                }
                            }
                        }
                        
                        // Checks if is into an aromatic ring
                        boolean isAromatic = false;
                        for (int k=0; k<nSK; k++) {
                            if ((k==j)||(k==i))
                                continue;
                            if (ConnMatrix[k][j] == 1.5) {
                                isAromatic = true;
                                break;
                            }
                        }
                        if (isAromatic)
                            nAromLinked++;
                    }
                }
                
                // counts H
                try {
                    nH = Mol.GetStructure().getAtom(i).getImplicitHydrogenCount();
                } catch (InvalidMoleculeException e) {
                    nH = 0;
                }
                
                // formal charge
                try {
                    Charge = Mol.GetStructure().getAtom(i).getFormalCharge();
                } catch (InvalidMoleculeException e) {
                    Charge = 0;
                }
                
                
                //// N_067: Al2-NH

                if ((Charge == 0) && (VD==2) && (nH == 1) && (nArom == 0))
                    if ((nAromLinked == 0) && (nSingle == 2) && (sOR == 0) && (sRdblX == 0)) 
                        N_067++;
                
                
                //// N_078: Ar-N=X / X-N=X

                if ((nTriple == 0) && (nElNeg >= 1) && (dElNeg >= 1) && (nO < 2))
                    N_078++;
                

                //// nRNNOx: Al2-N-N=O
                
                if ((nArom==0) && ((nH + sC)==2) && (sNdblO == 1))
                    nRNNOx++;
                        
                
                continue;
            }
            
            
            // Atom: O
            if (ConnMatrix[i][i] == 8) {
                
                
                continue;
            }
            
            
            // Atom: P
            if (ConnMatrix[i][i] == 15) {
                
                int sSorO=0, dSorO=0;
                
                for (int j=0; j<nSK; j++) {
                    if (j==i) 
                        continue;
                    if (ConnMatrix[i][j]>0) {
                        double Z = ConnMatrix[j][j];
                        double b = ConnMatrix[i][j];
                        // counts O or S neighbour
                        if ((Z==8)||(Z==16)) {
                            if (b==1)
                                sSorO++;
                            else if (b==2)
                                dSorO++;
                        }
                    }
                }
                
                if ((sSorO==3) && (dSorO==1))
                    nPO4++;
                
                continue;
            }
            
        }
        
    }

    /**
     * @return the nRNNOx
     */
    public int getnRNNOx() {
        return nRNNOx;
    }

    /**
     * @return the nPO4
     */
    public int getnPO4() {
        return nPO4;
    }

    /**
     * @return the N_067
     */
    public int getN_067() {
        return N_067;
    }

    /**
     * @return the N_078
     */
    public int getN_078() {
        return N_078;
    }
    
    
}
