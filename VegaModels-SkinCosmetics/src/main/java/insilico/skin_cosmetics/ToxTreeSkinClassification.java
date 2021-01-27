package insilico.skin_cosmetics;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.CustomQueryMatcher;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

import java.util.ArrayList;

/**
 *
 * @author Alberto
 */
public class ToxTreeSkinClassification {
    
    private final static String[][] SKIN_CLASSES = {
        // 0
        {"SNAr-Nucleophilic Aromatic Substitution", "identify a potential substrate for Nucleophilic aromatic substitution."}, 
        // 1
        {"Schiff Base Formation", "to be a potential substrate for Schiff base formation."}, 
        // 2
        {"Michael Acceptor", "to be a potential Michael reaction acceptor or a precursor."}, 
        // 3
        {"Acyl Transfer Agents", "identifies a potential acyl transfer agent or precursor."}, 
        // 4
        {"SN2-Nucleophilic Aliphatic Substitution", "to be a substrate for second-order nucleophilic substitution, or a precursor for it."}, 
        // 5
        {"No class found", "no alert for any sensitization class has been found."}, 
    };

    private final static String[][] SKIN_CLASSES_SMARTS = {
        // 0
        { 					
            "c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1",
            "c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])",
            "c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1",
            "c1([F,Cl,Br,I,$(N(=O)~O)])ncccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])",
            "c1([F,Cl,Br,I,$(N(=O)~O)])ncccn1",
            "c1([F,Cl,Br,I,$(N(=O)~O)])ncncc1",
            "c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])nc1",
            "c1nc([F,Cl,Br,I,$(N(=O)~O)])ncn1"
        },
        // 1
        {
            "[CH2]=O",
            "[CH2]N([CH3])[CH3]",
            "CC(C)=[CH][CH2][OH]",
            "CC=C(C)[CH2][OH]",
            "[CX4][CH]=[O,SX2]",
            "[a][CH]=O",
            "C(C)(C)=CC=[O,SX2]",
            "[C;!r5]([C;!r5])=[C;!r5](C)[C;!r5]=[O,SX2;!r5]",
            "[#6]C(=[O,SX2])C(=[O,SX2])[#6]",
            "[#6]C(=[O,SX2])[CX4]C(=[O,SX2])[#6]",
            "[#6][$([NX2]=O),$(N=C=O),$(OC#N),$(SC#N),$(N=C=S)]",
            "[CH2][NH2+0]"
        },
        // 2
        {
            "[CH2,CH]=[CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "C1=[C,N][$(S(=O)(=O)),$(C=[N,O]),$(S=O)][C,N]=C1c2ccccc2",
            "[CH2,CH]=C1C(=[O,SX2])**C1",
            "c1c([OH,NH2,NH])c([OH,NH2,NH,$(N=N),$(N(C)C)])ccc1",
            "c1c([OH,NH2,NH])ccc([OH,NH2,NH,$(N=N),$(N(C)C)])c1",
            "c1([OH])c(O[CH3])cccc1",
            "c1([OH])ccc(O[CH3])cc1",
            "c1c([OH])ccc(C=C[CH3])c1",
            "c1c([OH,NH2,NH])cc([OH,NH2,NH,$(N(C)C)])cc1",
            "C([F,Cl,Br,I])[CH2,CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "[$([CH]=[CH2,CH]),$(C(C)=[CH2,CH]),$(C#C);!$(C(C)=CC)][CH2][OH]",
            "[CH2]=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "c1([OH])ccc([CH2][CH]=[CH2])cc1",
            "c1(O[CH3])ccc([CH2][CH]=[CH2])cc1",
            "[OH]c1cccc2ccccc12",
            "F[CH2,CH]C(F)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "Cl[CH2,CH]C(Cl)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "Br[CH2,CH]C(Br)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "I[CH2,CH]C(I)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "aC=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "C#C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "[CH2,CH]=[CH]C=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "[CHR]=[CR][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]",
            "C=C([F,Cl,Br,I])[F,Cl,Br,I]",
            "c1c(=[O,NH2,NH])c(=[O,NH2,NH])ccc1",
            "c1c(=[O,NH2,NH])ccc(=[O,NH2,NH])c1"
        },
        // 3
        {
            "[!$(C=C);!$(C#C)]C(=[O,SX1,N])[F,Cl,Br,I]",
            "[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][a]",
            "[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N]C(=[O,SX1,N])",
            "[!$(C=C);!$(C#C)]C(=[O,S,N])[O,S][$(C=C),$(C=N),$(C#C),$(C#N)]",
            "[#6]1C(=[O,S])[O,N,S][#6]1",
            "[C,N,O,S,a]c1[n,o,s]c2ccccc2[n,o,s]1",
            "[#6]1[#6](=[N,O,S])[#7,#8,#16][#6][#7,#8,#16]1",
            "[C,O]=[#6]1[#7,#8,#16][#6](=[O,N,SX1])c2ccccc12",
            "[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][CX4,O,S][$(C=O),a,$(C=C),$(C#C),$(C=N),$(C#N)]"
        },
        // 4
        {
            "[CH,CH2,CH3;!$([CH2]CC=[O,S])][F,Cl,Br,I,$(OS(=O)(=O)[#6,#1]),$(OS(=O)(=O)O[#6,#1])]",
            "[#6]1[O,N,SX2][#6]1",
            "C1C(=[O,S])[O,S][CH2,CH]1",
            "[$(C=C),$(C#C),a][CH2,CH][O,S][$([CH]=O),$([CH]=S),$(C(C)=O),$(C(C)=S),a]",
            "[#6]1=,:[#6]C(=[O,SX1])N[SX2]1",
            "[#6][O,SX2,N][O,SX2,N][$([CH]=O),$([CH]=S),$(C(C)=O),$(C(C)=S),a]",
            "[CH2,CH3][NX3][NX2]=[O,S]",
            "c1ccc2cc3ccccc3cc2c1",
            "c1ccc2c(c1)ccc3ccccc23"
        },
        // 5
        { }
    };
    
    public class SkinClass {
        private final int Id;
        private final String Name;
        private final String Description;
        private final ArrayList<Pattern> Queries;
        
        public SkinClass(int ClassIndex) throws InitFailureException {
            this.Id = ClassIndex;
            this.Name = SKIN_CLASSES[ClassIndex][0];            
            this.Description = SKIN_CLASSES[ClassIndex][1];
            this.Queries = new ArrayList<>();
            for (String s : SKIN_CLASSES_SMARTS[ClassIndex])
                try {
                    Queries.add(SmartsPattern.create(s, DefaultChemObjectBuilder.getInstance()).setPrepare(false));
                } catch (Exception e) {
                    throw new InitFailureException("Unable to init SMARTS: " + s);
                }
        }
        
        public int getId() { return this.Id; }
        public String getName() { return this.Name; }
        public String getDescription() { return this.Description; }
        public ArrayList<Pattern> getQueries() { return this.Queries; }
    }
    
    
    
    public final static short SKIN_SNAR = 0;
    public final static short SKIN_SCHIFF = 1;
    public final static short SKIN_MICHAEL = 2;
    public final static short SKIN_ACYL = 3;
    public final static short SKIN_SN2 = 4;
    public final static short SKIN_NO_CLASS = 5;
    
    private final SkinClass SNAr;
    private final SkinClass SchiffBase;
    private final SkinClass MichaelAcceptor;
    private final SkinClass AcylTransfer;
    private final SkinClass SN2;
    private final SkinClass NoClass;

    
    public ToxTreeSkinClassification() throws InitFailureException {
        // Init class rules
        SNAr = new SkinClass(SKIN_SNAR);
        SchiffBase = new SkinClass(SKIN_SCHIFF);
        MichaelAcceptor = new SkinClass(SKIN_MICHAEL);
        AcylTransfer = new SkinClass(SKIN_ACYL);
        SN2 = new SkinClass(SKIN_SN2);
        NoClass = new SkinClass(SKIN_NO_CLASS);
    }
    
    
    public SkinClass Calculate(InsilicoMolecule mol) throws GenericFailureException {
        
        CustomQueryMatcher Matcher;
        try {
            Matcher = new CustomQueryMatcher(mol);
        } catch (CDKException e) {
            throw new GenericFailureException("Unable to init SMARTS matcher");
        }
        
        if (Matches(mol, SNAr.getQueries()))
            return SNAr;
        if (Matches(mol, SchiffBase.getQueries()))
            return SchiffBase;
        if (Matches(mol, MichaelAcceptor.getQueries()))
            return MichaelAcceptor;
        if (Matches(mol, AcylTransfer.getQueries()))
            return AcylTransfer;
        if (Matches(mol, SN2.getQueries()))
            return SN2;
        
        return NoClass;
    }
    
    
    private boolean Matches(InsilicoMolecule mol, ArrayList<Pattern> Queries )
            throws GenericFailureException {
        try {
            for (Pattern q : Queries)
                if (q.matches(mol.GetStructure()))
                    return true;
        } catch (InvalidMoleculeException e) {
            throw new GenericFailureException("Error while performing matching - " + e.getMessage());
        }
        return false;
    }

}
