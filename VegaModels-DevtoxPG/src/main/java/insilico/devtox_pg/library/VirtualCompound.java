package insilico.devtox_pg.library;

import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import java.io.Serializable;
import java.util.BitSet;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

/**
 *
 * @author Alberto
 */
public class VirtualCompound implements Serializable {

    private String SMILES;
    private String SMARTS;
    private int nAtoms;
    private Pattern Query;
    private BitSet FP;
    
    
    public VirtualCompound(String MolSMILES, boolean KeepQuery) throws CDKException {
        
        InsilicoMolecule curMol = SmilesMolecule.Convert(MolSMILES);
        
        this.SMILES = MolSMILES;
        try {
            this.nAtoms = curMol.GetStructure().getAtomCount();
        } catch (Throwable ex) {
            throw new CDKException("Unable to get structure - " + ex.getMessage());
        }

        // pass structure through custom smiles/smarts writer
        // so to avoid matching on atoms which don't belong to
        // any ring, this is set as the SMARTS for the compound
        try {
            CustomSmilesWriterWithRings sg = new CustomSmilesWriterWithRings();
            sg.setUseAromaticityFlag(true);
            this.SMARTS = sg.createSMILES(curMol.GetStructure());
        } catch (Throwable ex) {
            throw new CDKException("Unable to build smarts - " + ex.getMessage());
        }

        if (KeepQuery) {
            try {
                this.Query = SmartsPattern.create(this.SMARTS, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            } catch (Throwable ex) {
                throw new CDKException("Unable to parse as SMARTS the string: " + MolSMILES);
            }
        } else
            this.Query = null;
        
        // Calculate FP for the query
        try {
            AllRingsFinder ARF = new AllRingsFinder();
            ARF.setTimeout(10000);
            Fingerprinter fp = new Fingerprinter();
            FP = fp.getFingerprint(curMol.GetStructure());
        } catch (Throwable ex) {
            throw new CDKException("Unable to create FP: " + ex.getMessage());
        }
            
    }
    
    /**
     * @return the SMILES
     */
    public String getSMILES() {
        return SMILES;
    }

    /**
     * @return the SMARTS
     */
    public String getSMARTS() {
        return SMARTS;
    }

    /**
     * @return the nAtoms
     */
    public int getnAtoms() {
        return nAtoms;
    }

    /**
     * @return the Query
     */
    public Pattern getQuery() {
        return Query;
    }

    /**
     * @return the QueryFP
     */
    public BitSet getFP() {
        return FP;
    }       
    
}
