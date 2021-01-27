package insilico.fathead_epa.descriptors;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.CustomQueryMatcher;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

/**
 *
 * @author amanganaro
 */
public class LC50FMDescriptors {

    private int n_aNO2;
    private int n_CHO;
    private int n_POH;
    private int n_NNO;
    private int n_NCON;
    private int n_CN;
    private int n_Cl;
    private int n_aCHO;
    private int n_sssS;
    
    Pattern[] q;
    
    
    public LC50FMDescriptors() throws InitFailureException {

        q = new Pattern[9];
        try {
            q[0] = SmartsPattern.create("aN(~O)~O").setPrepare(false);
            q[1] = SmartsPattern.create("[!a][C;D2]=O").setPrepare(false);
            q[2] = SmartsPattern.create("P[O;D1]").setPrepare(false);
            q[3] = SmartsPattern.create("[N;D3](-*)(-*)[N;D2]=O").setPrepare(false);
            q[4] = SmartsPattern.create("NC(=O)N").setPrepare(false);
            q[5] = SmartsPattern.create("[!a]C#N").setPrepare(false);
            q[6] = SmartsPattern.create("C=C-Cl").setPrepare(false);
            q[7] = SmartsPattern.create("a[C;D2]=O").setPrepare(false);
            q[8] = SmartsPattern.create("[S;D3](-*)(-*)-*").setPrepare(false);
        } catch (Exception ex) {
            throw new InitFailureException("unable to init SMARTS parsers");
        }
        
    }
    
    
    public void Calculate(InsilicoMolecule mol) throws GenericFailureException {
        
        try {

            CustomQueryMatcher Matcher = new CustomQueryMatcher(mol);
            
            n_aNO2 = 0;
            if (q[0].matches(mol.GetStructure())) {
                n_aNO2 = q[0].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_CHO = 0;
            if (q[1].matches(mol.GetStructure())) {
                n_CHO = q[1].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_POH = 0;
            if (q[2].matches(mol.GetStructure())) {
                n_POH = q[2].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_NNO = 0;
            if (q[3].matches(mol.GetStructure())) {
                n_NNO = q[3].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_NCON = 0;
            if (q[4].matches(mol.GetStructure())) {
                n_NCON = q[4].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_CN = 0;
            if (q[5].matches(mol.GetStructure())) {
                n_CN = q[5].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_Cl = 0;
            if (q[6].matches(mol.GetStructure())) {
                n_Cl = q[6].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_aCHO = 0;
            if (q[7].matches(mol.GetStructure())) {
                n_aCHO = q[7].matchAll(mol.GetStructure()).countUnique();
            }
            
            n_sssS = 0;
            if (q[8].matches(mol.GetStructure())) {
                n_sssS = q[8].matchAll(mol.GetStructure()).countUnique();
            }
            
            
        } catch (Exception e) {
            throw new GenericFailureException("unable to perform matching - " + e.getMessage());
        }
                
        
    }

    /**
     * @return the n_aNO2
     */
    public int getN_aNO2() {
        return n_aNO2;
    }

    /**
     * @return the n_CHO
     */
    public int getN_CHO() {
        return n_CHO;
    }

    /**
     * @return the n_POH
     */
    public int getN_POH() {
        return n_POH;
    }

    /**
     * @return the n_NNO
     */
    public int getN_NNO() {
        return n_NNO;
    }

    /**
     * @return the n_NCON
     */
    public int getN_NCON() {
        return n_NCON;
    }

    /**
     * @return the n_CN
     */
    public int getN_CN() {
        return n_CN;
    }

    /**
     * @return the n_Cl
     */
    public int getN_Cl() {
        return n_Cl;
    }

    /**
     * @return the n_aCHO
     */
    public int getN_aCHO() {
        return n_aCHO;
    }

    /**
     * @return the n_sssS
     */
    public int getN_sssS() {
        return n_sssS;
    }
     
}
