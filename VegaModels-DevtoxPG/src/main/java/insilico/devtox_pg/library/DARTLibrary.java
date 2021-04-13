package insilico.devtox_pg.library;

import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.similarity.Similarity;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.net.URL;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

/**
 *
 * @author Alberto
 */
public class DARTLibrary {
    
    private final static boolean USE_SCAFFOLD = false;

    private final static Object[][] CategoryNames = {
        {1, "Inorganics and derivatives: metals, metallic derivatives, organophosphorus" +
            "and organosiloxane compounds"},
        {2, "Estrogen receptor (ER) and androgen receptor (AR) binding compounds"},
        {3, "Retinoic acid receptor (RAR), aryl hydrocarbon receptor (AhR) binders and " +
            "prostaglandin receptor agonists"},
        {4, "Nicotinic acetyl choline receptor (nAChR) binders and acetyl cholineesterase " +
            "(AChE) inhibitors"},
        {5, "Ion channel opener/inhibitor, beta-adrenergic inhibitors, ACE/ARA inhibitors " +
            "and Shh signaling interference/Cholesterol synthesis inhibition"},
        {6, "Opioid receptor binders and tubulin receptor interactors"},
        {7, "Nucleotide and nucleobase derivatives"},
        {8, "Aromatic compounds with alkyl, multi-halogen and nitro groups"},
        {9, "Aromatic compounds contain alkyl chain with alcohol, aldehyde, acid " +
            "functional groups; poly-Cl aryloxy derived acids and esters"},
        {10, "Aromatic compounds with sulfonamide and urea moieties, phenytoins"},
        {11, "Aromatic compounds (non-fused ring system) with aliphatic " +
            "amine moieties"},
        {12, "Aromatic diamine, their diazo moieties, and aromatic triazene " +
            "derivatives"},
        {13, "Imidazole, nitro imidazoles derivatives, nitro-furfurylideneamino and triazole " +
            "derivatives"},
        {14, "Aromatic ring fused cyclic-, heterocyclic derivatives"},
        {15, "Miscellaneous aromatic chemicals and antibiotics"},
        {16, "Non-aromatic cyclic hydrocarbon ring, heterocyclic ring contain only oxgen" +
            "atom and multi Cl single/fused cyclic hydrocarbons"},
        {17, "Heterocyclic, cyclic compounds contain nitrogen, oxygen/sulfur atoms"},
        {18, "Miscellaneous cyclic chemicals"},
        {19, "Alkyl carbamodi-thioic acids, alkyl sulfonates and perfluorinated compounds " +
            "(PFCs)"},
        {20, "Miscellaneous non-cyclic chemicals"},
        {21, "Vinyl amides, aldehydes, esters and alkyl amides (<C4), N-substituted " +
            "amides, ureas, carbonates, guanidine and carbamates"},
        {22, "Alpha-substituted carboxylic acids, esters and di-acid derived esters"},
        {23, "Small- (C1-C4) halo-, multihalo-alkanes, alkyl ether/alkenes and halogenated " +
            "acetonitriles as well as N, or S related mustards"},
        {24, "Di/multi-OH, NH2, substituted amine, SH (=S), OR, OAc substituted (at each " +
            "terminal carbon) C1 to C5 hydrocarbon chain or repeating C2 units as well as metal " +
            "chelators"},
        {25, "C1 to C4 non-branched/<C9 beta-alkyl (<C5) substituted alcohols and <C4 " +
            "alkyl, vinyl nitriles"},
        {26, "Inactive compounds"}
    };    
    
    // queries for rule 01
    private final Pattern s1_a_1, s1_a_2_a, s1_a_2_b, s1_a_3;
    private final Pattern s1_b;
    private final Pattern s1_c_1, s1_c_2;
    private final Pattern s_hetero;
    private final Pattern s8a_toluene;
    private final Pattern s2b_3_3a, s2b_3_3b, s2b_3_3c;
    private final InsilicoMolecule[] s20;
    
    private static final String LIB_ROOT = "/dart_library/";
    private Pattern CurMatcher;
    private final AllRingsFinder ARF;
    private final Fingerprinter Fprinter;
    
    
    public class DARTResult {
        public int Index;
        public String SubIndex;
        public boolean isScaffold;
        public String Structure;
        
        public DARTResult(int Index, String SubIndex, boolean isScaffold, String Structure) {
            this.Index = Index;
            this.SubIndex = SubIndex;
            this.isScaffold = isScaffold;
            this.Structure = Structure;
        }
    }
    
    
    public DARTLibrary() throws InitFailureException {
        ARF = new AllRingsFinder();
//        ARF.setTimeout(10000);
        Fprinter = new Fingerprinter();   
        
        try {
            //// 1a - Metals and metallic derivative
            // As, B, Mn, Cr, Zn, Te
            s1_a_1 = SmartsPattern.create("[#5,#33,#25,#24,#30,#52]").setPrepare(false);

            // Al, Cd, Cu, Zn, Mn, Ni, Pb chlorides
            s1_a_2_a = SmartsPattern.create("[#13,#48,#29,#30,#25,#28,#82]").setPrepare(false);
            s1_a_2_b = SmartsPattern.create("[Cl-]").setPrepare(false);

            // Pb, Hg, Sn derivatives
            s1_a_3 = SmartsPattern.create("[#82,#80,#50]").setPrepare(false);
            
            //// 1b - Organophosphorus compounds
            s1_b = SmartsPattern.create("P=[O,S]").setPrepare(false);
            
            //// 1c - Organosiloxane compounds
            // Organosiloxanes
            s1_c_1 = SmartsPattern.create("[O;D2][Si]([C;D1])[O;D2]").setPrepare(false);
            
            // Phenyl-siloxanes
            s1_c_2 = SmartsPattern.create("c1ccccc1[Si]([C;D1])[O;D2][Si]([C;D1])").setPrepare(false);

            //// 2b-3-3 a, b and c
            s2b_3_3a = SmartsPattern.create("O=C([O;!$(O(C([c])=O)C[!C]);!$(O(C([c])=O)CC[!C]);!$(O(C([c])=O)CCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCCC[!C]);!$(O(C([c])=O)CCCCCC[!C])])c1ccc([OH])cc1").setPrepare(false);
            s2b_3_3b = SmartsPattern.create("O=C([O;!$(O(C([c])=O)C[!C]);!$(O(C([c])=O)CC[!C]);!$(O(C([c])=O)CCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCCC[!C]);!$(O(C([c])=O)CCCCCC[!C])])c1ccccc1(C(=O)[O;!$(O(C([c])=O)C[!C]);!$(O(C([c])=O)CC[!C]);!$(O(C([c])=O)CCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCCC[!C]);!$(O(C([c])=O)CCCCCC[!C])])").setPrepare(false);
            s2b_3_3c = SmartsPattern.create("O=C([O;!$(O(C([c])=O)C[!C]);!$(O(C([c])=O)CC[!C]);!$(O(C([c])=O)CCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCC[!C]);!$(O(C([c])=O)CCCCC[!C]);!$(O(C([c])=O)CCCCCC[!C])])c1ccccc1([OH])").setPrepare(false);
            
            //// 8a - smarts for toluene
            s8a_toluene = SmartsPattern.create("[CH3]c1ccccc1").setPrepare(false);
            
            //// 20 - all sub classes
            s20 = new InsilicoMolecule[8];
            s20[0] = SmilesMolecule.Convert("O=C(OCN=[N+]([O-])C)C");
            s20[1] = SmilesMolecule.Convert("CCCCCC");
            s20[2] = SmilesMolecule.Convert("O=C(C)CCCC");
            s20[3] = SmilesMolecule.Convert("O=C(C)CCC(=O)C");
            s20[4] = SmilesMolecule.Convert("O=C(CCl)C(Cl)(Cl)Cl");
            s20[5] = SmilesMolecule.Convert("O=C(C(Cl)Cl)C(Cl)Cl");
            s20[6] = SmilesMolecule.Convert("O=C(C(F)(F)F)C(F)(F)F");
            s20[7] = SmilesMolecule.Convert("O=C" );
           
            //// check for heteroatoms (for PAH and toluene)
            s_hetero = SmartsPattern.create("[!#6]").setPrepare(false);

            
        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS");
        }           
    }
    
    
    public DARTResult Check(InsilicoMolecule mol) throws IOException, ClassNotFoundException, InvalidMoleculeException, CDKException {
        

        //// check category 01
        
        try {

            if ( (s1_a_1.matches(mol.GetStructure())) || ( (s1_a_2_a.matches(mol.GetStructure()))) && ((s1_a_2_b).matches(mol.GetStructure())) || ((s1_a_3).matches(mol.GetStructure())) )
                return(new DARTResult(1, "a", false, "Metals and metallic derivatives"));            
            if ( ((s1_b).matches(mol.GetStructure()) ))
                return(new DARTResult(1, "b", false, "Organophosphorus compounds"));            
            if ( ((s1_c_1).matches(mol.GetStructure()) || ((s1_c_2).matches(mol.GetStructure()) )))
                return(new DARTResult(1, "c", false, "Organosiloxane compounds"));            
        } catch (Exception e) {
            throw new CDKException("unable to check category 01");
        }        
        
        
        //// check (automatic) all remaining categories
        
        DARTResult ResCategory = null;
        DARTResult ResScaffold = null;
        int ResScaffoldAtoms = 0;
        int nAt = mol.GetStructure().getAtomCount();
        CurMatcher = SmartsPattern.create(mol.GetSMILES()).setPrepare(false);

        // calculate FP for the molecule
        BitSet CurFP;
        try {
            CurFP = Fprinter.getFingerprint(mol.GetStructure());
        } catch (Throwable ex) {
            throw new CDKException("unable to calculate FP - " + ex.getMessage());
        }
        
        // libraries of compounds 2-26
        for (int lib=2; lib<27; lib++) {
            
            // Hardcoded exception for rule 20
            // Force exact match with compounds
            if ( (lib == 20) ) {
                for (InsilicoMolecule vc : s20)
                    if (IsExactMatching(mol, vc)) {
                        ResCategory = new DARTResult(20, "", false, vc.GetSMILES());
                        break;                                                            
                    }
                // no matches - skip lib
                continue;
            } 
            
            // de-serialize current library
            URL LibURL = getClass().getResource(LIB_ROOT + lib + ".dat");            
            ObjectInputStream in = new ObjectInputStream(LibURL.openStream());
            VirtualCompoundCategory VCC = (VirtualCompoundCategory) in.readObject();
            in.close();    
            
            // check library
            for (VirtualCompoundSubCategory SVCC : VCC.getSubCategory()) {
                
                boolean ProceedToMolecules = true;
                
                // in current version, scaffolds are not checked
                if (USE_SCAFFOLD) {
                    // if scaffolds are available, check SMARTS matching on them
                    // otherwise does not proceed to single virtual compounds
                    if (SVCC.getScaffoldsSize() > 0) {
                        ProceedToMolecules = false;
                        for (int i=0; i<SVCC.getScaffoldsSize(); i++) {

                            Pattern pattern = SmartsPattern.create(SVCC.getScaffold(i).getQuery().toString()).setPrepare(false);
                            if (pattern.matches(mol.GetStructure())) {
                                ProceedToMolecules = true;
                                if (SVCC.getScaffold(i).getnAtoms() > ResScaffoldAtoms) {
                                    ResScaffold = new DARTResult(SVCC.getIndex(), SVCC.getSubIndex(), true, SVCC.getScaffold(i).getSMILES());
                                    ResScaffoldAtoms = SVCC.getScaffold(i).getnAtoms();
                                }
                                break;
                            }
                        }
                    }
                }                

                // check single virtual compounds
                // proceed to matching only if the possible match will represent
                // more than 60% (no. of atoms) of the target molecule,
                // then if the FP subset is matching
                if (ProceedToMolecules) {

                    // Some hardcoded exception for specific categories
                    // Force exact match with toluene before proceeding
                    
                    // 2b-3-3a
                    if ( (SVCC.getIndex()==2) && (SVCC.getSubIndex().equalsIgnoreCase("b-3-3a") )) 
//                        if (!Matcher.matches(s2b_3_3a)) continue;
                        if (!(s2b_3_3a).matches(mol.GetStructure())) continue;

                    // 2b-3-3
                    if ( (SVCC.getIndex()==2) && (SVCC.getSubIndex().equalsIgnoreCase("b-3-3b") ))
                        if (!(s2b_3_3b).matches(mol.GetStructure())) continue;

                    // 2b-3-3c
                    if ( (SVCC.getIndex()==2) && (SVCC.getSubIndex().equalsIgnoreCase("b-3-3c") ))
                        if (!(s2b_3_3c).matches(mol.GetStructure())) continue;

                    // 3b-3 (skip PAH if heteroatoms found)
                    if ( (SVCC.getIndex()==3) && (SVCC.getSubIndex().equalsIgnoreCase("b-3") ))
                        if (!(s_hetero).matches(mol.GetStructure())) continue;
                        
                    // 8a
                    if ( (SVCC.getIndex()==8) && (SVCC.getSubIndex().equalsIgnoreCase("a") )) {
                        if (!(s8a_toluene).matches(mol.GetStructure())) continue;
                        if ((s_hetero).matches(mol.GetStructure())) continue;
                    }

                    if (SVCC.getMoleculesSize() > 0) {
                        for (int i=0; i<SVCC.getMoleculesSize(); i++) {
                            
                            double CurOverlap = (double)SVCC.getMolecule(i).getnAtoms() / (double)nAt;
                            if (CurOverlap > 0.6) {

                                boolean isSubstr = FingerprinterTool.isSubset(CurFP, SVCC.getMolecule(i).getFP());
                                if (isSubstr) {
                                    Pattern pattern = SmartsPattern.create(SVCC.getMolecule(i).getSMARTS()).setPrepare(false);
                                    if (pattern.matches(mol.GetStructure())) {
                                        ResCategory = new DARTResult(SVCC.getIndex(), SVCC.getSubIndex(), false, SVCC.getMolecule(i).getSMILES());
                                        break;                                    
                                    }
                                }
                            }
                            
                        }
                    }
                }
                
                if (ResCategory != null)
                    return ResCategory;
            }
            
            // free resource
            VCC = null;
        }
        
        if (ResScaffold != null ) {
            return (ResScaffold);
        } else
            return null;
    }    
    

//    boolean IsMatching(String VirtualCompoundSMARTS)  {
//
//        try {
//            QueryAtomContainer CurQuery = SMARTSParser.parse(VirtualCompoundSMARTS);
//            if (CurMatcher.matches(CurQuery))
//                return true;
//        } catch (Throwable ex) {
//            System.out.println("Error on " + VirtualCompoundSMARTS);
//        }
//        return false;
//    }
    
    boolean IsMatching(QueryAtomContainer VirtualCompoundQuery)  {

        try {
            if (CurMatcher.matches(VirtualCompoundQuery))
                return true;
        } catch (Throwable ex) {
            //System.out.println("Error on " + VirtualCompoundSMILES);
        }
        return false;
    }
    
    boolean IsExactMatching(InsilicoMolecule Mol, InsilicoMolecule VirtualCompound)  {
        
        try {
            if (Similarity.CheckIsomorphism(Mol.GetStructure(),VirtualCompound.GetStructure()))
                return true;
        } catch (Throwable ex) {
            //System.out.println("Error on " + VirtualCompoundSMILES);
        }
        return false;
    }
    
    
    /**
     * Parse all SMARTS, only of virtual compounds, in all 2-26 libraries.
     * 
     * @throws IOException
     * @throws ClassNotFoundException 
     */
    public void CheckVirtualCompoundsSMARTS() throws IOException, ClassNotFoundException {
        
        // libraries of compounds 2-26
        for (int lib=2; lib<27; lib++) {
            
            // de-serialize current library
            URL LibURL = getClass().getResource(LIB_ROOT + lib + ".dat");            
            ObjectInputStream in = new ObjectInputStream(LibURL.openStream());
            VirtualCompoundCategory VCC = (VirtualCompoundCategory)in.readObject();
            in.close();    
            
            // check library
            for (VirtualCompoundSubCategory SVCC : VCC.getSubCategory()) {
                
                if (SVCC.getMoleculesSize() > 0) {
                    for (int i=0; i<SVCC.getMoleculesSize(); i++) {
                        String smi = SVCC.getMolecule(i).getSMILES();
                        String sma = SVCC.getMolecule(i).getSMARTS();
                        String id = SVCC.getIndex() + SVCC.getSubIndex();
                        try {
                            Pattern CurQuery = SmartsPattern.create(sma).setPrepare(false);
                        } catch (Throwable ex) {
                            System.out.println(id + "\t" + "unable to parse smarts:\t" + sma + "\t" + smi);
                        }
                            
                    }
                }
            }
        }
    }        
    
    
    final static public String getCategoryName(int Index) {
        
        for (int i=0; i<CategoryNames.length; i++) 
            if ( (int)CategoryNames[i][0] == Index )
                return (String)CategoryNames[i][1];        
        
        return "-";
    }
}
