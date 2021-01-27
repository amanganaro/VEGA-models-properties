package insilico.bcf_arnotgobas.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.FunctionalGroups;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.ringsearch.RingPartitioner;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.List;

/**
 * kM Biotransformation Rate in Fish - factor for the main formula
 * Calculates KM as implemented in EpiSuite.
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class KmFactor extends DescriptorBlock {

    private final static long serialVersionUID = 1L;
    private final static String BlockName = "kM Biotransformation Rate in Fish";

    private final static String[][] f_31 = {
            {"unsubstituted 5-membered ring", "[a;!$([a][A])]1[a;!$([a][A])][a;!$([a][A])][a;!$([a][A])][a;!$([a][A])]1"},
            {"unsubstituted 6-membered ring", "[a;!$([a][A])]1[a;!$([a][A])][a;!$([a][A])][a;!$([a][A])][a;!$([a][A])][a;!$([a][A])]1"},
    };

    private final static String Benzene = "c1ccccc1";

    // Array with: index, name, coefficient, SMARTS (if needed, null otherwise)
    private final static Object[][] KmFragments = {
            {1, "Nitroso [-N-N=O]", -0.42851048, "O=[N;!R]-[N;!R;$([N;D3]([C,a])[C,a]),$([N;D2][C,a]),$([N;D1])]"},
            {2, "Linear C4 terminal chain [CCC-CH3]", 0.03412373, "[C;D2][C;D2][C;D2][CH3;D1]"},
            {3, "Aliphatic alcohol [-OH]", -0.06155701, null},
            {4, "Aromatic alcohol [-OH]", -0.47273947, null},
            {5, "Aliphatic acid [-C(=O)-OH]", 0.38030117, null},
            {6, "Aldehyde [-CHO]", 0.24648749, null},
            {7, "Ester [-C(=O)-O-C]", -0.76052851, "O=C[O;D2][C,c]"},
            {8, "Amide [-C(=O)-N or C(=S)-N]", -0.59521049, "[O,S]=[C;!R;!$([C;D3](N)S);!$([C;D3](N)P);!$([C;D3](N)O);!$([C;D3](N)N)]-N"},
            {9, "Triazine ring (symmetric)", -0.01226285, null},
            {10, "Aliphatic chloride [-CL]", 0.36076089, "[Cl][A]"},
            {11, "Aromatic chloride [-CL]", 0.37784643, "[Cl][a]"},
            {12, "Aliphatic bromide [-Br]", 0.27340813, "[Br][A]"},
            {13, "Aromatic bromide [-Br]", 0.39635369, "[Br][a]"},
            {14, "Aromatic iodide [-I]", 0.21395406, "[I][a]"},
            {15, "Carbon with 4 single bonds & no hydrogens", -0.29842827, "[C;D4;!$(C(F)(F)F)]"},
            {16, "Aromatic nitro [-NO2]", -0.02177166, null},
            {17, "Aliphatic amine [-NH2 or -NH-]", 0.40673985, null},
            {18, "Aromatic amine [-NH2 or -NH-] ", -0.28895783, null},
            {19, "Cyanide / Nitriles [-C#N]", 0.1542211, null},
            {20, "Sulfonic acid / salt -> aromatic attach", 0.0247512, "[S;D4](=[S,O])(=[S,O])([S,O;D1])[a;!s;!o]"},
            {21, "Polyaromatic hydrocarbon (4 or more rings)", 0.46164904, "[C,a]"},
            {22, "Pyridine ring", -0.90212009, null},
            {23, "Aromatic ether [-O-aromatic carbon]", -0.0693912, "[c][O;D2][C,c]"},
            {24, "Aliphatic ether [C-O-C]", -0.02324019, "C-O-C"},
            {25, "Ketone [-C-C(=O)-C-]", -0.1800634, "[C,c]C(=O)[C,c]"},
            {26, "Tertiary amine", -0.78292477, null},
            {27, "Phosphate ester (P=O type)", -0.60313654, "[$(P([O;D2])=O)]"},
            {28, "Alkyl substituent on aromatic ring", 0.17805958, "[a][C;!$(C~N);!$(C~O);!$(C~P);!$(C~S);!$(C~[Cl]);!$(C~[Br]);!$(C~[I]);!$(C~[F]);!$(C([C;D1])([C;D1])[C;D1])]"},
            {29, "Carbamate or Thiocarbamate", -0.77660157, "C(=[O,S])([$([O;D1]),$([S;D1]),$([O;D2](C)C),$([S;D2](C)C)])[$([N;D1]),$([N;D2](C)[C,c]),$([N;D3](C)([C,c])[C,c])]"},
            {30, "Trifluoromethyl group [-CF3]", -0.18814672, "C(F)(F)F"},
            {31, "Unsubstituted aromatic (3 or less rings)", 0.10024087, null},
            {32, "Unsubstituted phenyl group (C6H5-)", -0.60319946, "[*]-[c;D3]1[c;D2][c;D2][c;D2][c;D2][c;D2]1"},
            {33, "Phosphate ester (P=S type)", 0.19779672, "[$(P([O;D2])=S)]"},
            {34, "Urea [N-C(=O)-N]", -1.25487927, null},
            {35, "Furan or Thiofuran", -0.35707841, null},
            {36, "Triazole Ring ", 0.32253333, null},
            {37, "Fluorine  [-F]", 0.27586509, "[F]"},
            {38, "Aromatic-CH3", -0.08716123, "[CH3;D1;$(C[a])]"},
            {39, "Aromatic-CH2", -0.33650743, "[CH2;D2;$(C[a])]"},
            {40, "Aromatic-CH", -0.46293351, "[CH1;D3;$(C[a])]"},
            {41, "Aromatic-H", 0.26637806, "[a;H1]"},
            {42, "Methyl [-CH3]", 0.24510529, "[A][CH3;D1]"},
            {43, "-CH2- [linear]", 0.02418707, "[CH2;D2;!R;!R2;!R3;!$(C[a])]"},
            {44, "-CH- [linear]", -0.19123158, "[CH;D3;!R;!R2;!R3;!$(C[a])]"},
            {45, "-CH2- [cyclic]", 0.09625069, "[CH2;D2;R,R2,R3]"},
            {46, "-CH - [cyclic]", 0.01260466, "[CH;D3;R,R2,R3]"},
            {47, "-C=CH [alkenyl hydrogen]", 0.09884729, "C=[CH]"},
            {48, "Thiazole Ring", -0.44077474, null},
            {49, "o-Chloro / Mono-aromatic ether", 0.78281744, "C[O;D2][$(aa[Cl])]"},
            {50, "Number of fused acyclic rings", 0.647663, null},
            {51, "Number of fused 6-carbon aromatic rings", -0.5778540, "[$(c1ccc2ccccc2c1)]1ccccc1"},
            {52, "Four or more fused aromatic rings", -0.512825, null},
            {53, "Four or more fused cyclic rings ", -1.727859, null},
            {54, "Benzene", -0.427728, "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1"},
            {55, "Naphthalene", 0.431984, "[c;R1]1[c;R1][c;R1]c2[c;R1][c;R1][c;R1][c;R1]c2[c;R1]1"},
            {56, "Indane", -0.454496, "c1ccc2c(c1)CCC2"},
            {57, "Biphenyl", -0.531935, "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]2"},
    };

//    private final QueryAtomContainer[] Queries; // one for each pre-fetched SMARTS
//    private final QueryAtomContainer[] Queries_f_31; // for fragment 31
//    private final QueryAtomContainer Query_benzene; // just benzene ring

    private Pattern[] Queries; // one for each pre-fetched SMARTS
    private Pattern[] Queries_f_31; // for fragment 31
    private Pattern Query_benzene; // just benzene ring

    // Private class to perform ring analysis for some fragments
    private class RingAnalyzer {

        public int nFusedBenzene;
        public int nFusedAliphatic;
        public int nRings;
        public boolean FourFusedAromatic;
        public boolean FourFusedAliphatic;

        public RingAnalyzer(InsilicoMolecule mol) throws Exception {

            nFusedBenzene = 0;
            nFusedAliphatic = 0;
            FourFusedAromatic = false;
            FourFusedAliphatic = false;

            // Detects rings
            IRingSet sssrings;
            sssrings = mol.GetSSSR();
            nRings = sssrings.getAtomContainerCount();

            if (nRings > 1) {

                // partition ring into fused rings sets
                List<?> ringsets = RingPartitioner.partitionRings(sssrings);

                for (Object o : ringsets) {

                    IRingSet ringset = (IRingSet) o;

                    if (ringset.getAtomContainerCount() < 2)
                        continue;

                    int aliph_ring_count = 0;
                    int bnz_ring_count = 0;

                    for (int j = 0; j < ringset.getAtomContainerCount(); j++) {

                        boolean isArom = true;
                        boolean isBnz = true;
                        if (ringset.getAtomContainer(j).getAtomCount() != 6)
                            isBnz = false;
                        for (IAtom a : ringset.getAtomContainer(j).atoms()) {
                            if (!a.getFlag(CDKConstants.ISAROMATIC))
                                isArom = false;
                            if (a.getAtomicNumber() != 6)
                                isBnz = false;
                        }

                        if ((isArom) && (isBnz)) bnz_ring_count++;
                        if (!isArom) aliph_ring_count++;
                    }


                    boolean FourArom = false;
                    boolean FourAliph = false;

                    if (ringset.getAtomContainerCount() >= 4) {
//                        if (bnz_ring_count >= 4)
                        if (bnz_ring_count >= aliph_ring_count)
                            FourArom = true;
                        else
                            FourAliph = true;
                    }

                    if (FourArom) FourFusedAromatic = true;
                    if (FourAliph) FourFusedAliphatic = true;

                    nFusedBenzene += bnz_ring_count;

                    if (ringset.getAtomContainerCount() >= 4) {
                        if (bnz_ring_count >= 4) {
                            if (aliph_ring_count > 1)
                                nFusedAliphatic++;
                        } else if (aliph_ring_count >= 4) {
                            nFusedAliphatic += aliph_ring_count;
                        } else {
                            if (aliph_ring_count > 1)
                                nFusedAliphatic += aliph_ring_count;
                            //nFusedAliphatic++;
                        }
                    } else if (aliph_ring_count > 1)
                        nFusedAliphatic++;

                }
            }
        }
    }


    /**
     * Constructor. This should not be used, no weight is specified. The
     * overloaded constructors should be used instead.
     */
    public KmFactor() {
        super();
        this.Name = KmFactor.BlockName;

        // Init needed SMARTS
        Queries = new Pattern[KmFragments.length];
        for (int i = 0; i < KmFragments.length; i++) {
            if (KmFragments[i][3] != null)
                try {
                    Queries[i] = SmartsPattern.create((String) KmFragments[i][3]).setPrepare(false);
                } catch (Exception e) {
                    Queries[i] = null;
                }
        }

        Queries_f_31 = new Pattern[f_31.length];
        for (int i = 0; i < f_31.length; i++) {
            try {
                Queries_f_31[i] = SmartsPattern.create(f_31[i][1]).setPrepare(false);
            } catch (Exception e) {
                Queries_f_31[i] = null;
            }
        }


        try {
            Query_benzene = SmartsPattern.create(Benzene).setPrepare(false);
        } catch (Exception e) {

        }

    }


    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        Add("kM", "kM Biotransformation Rate in Fish");

        // Create all single fragments
        for (int i = 0; i < KmFragments.length; i++)
            Add("kM_f_" + (i + 1), "kM fragment " + (i + 1) + ": " + KmFragments[i][1]);

        SetAllValues(Descriptor.MISSING_VALUE);
    }


    /**
     * Calculate descriptors for the given molecule.
     *
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {
        this.Calculate(mol, null);
    }


    /**
     * Overload of Calculate(), to be used if needed descriptors are passed (to
     * fasten calculation)
     *
     * @param mol         molecule to be calculated
     * @param FunctGroups Functional groups descriptors
     */
    public void Calculate(InsilicoMolecule mol, DescriptorBlock FunctGroups) {

        // Generate/clears descriptors
        GenerateDescriptors();

        IAtomContainer molStructure;
        try {
            molStructure = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        // Retrieves or calculates Functional groups
        DescriptorBlock CurFG;
        if (FunctGroups != null)
            CurFG = FunctGroups;
        else {
            CurFG = new FunctionalGroups();
            CurFG.Calculate(mol);
        }

        // Calculates or derive all fragments

        try {

            // Matcher tool from old insilico
//            CustomQueryMatcher Matcher;
//            Matcher = new CustomQueryMatcher(mol);

            // Ring analysis
            RingAnalyzer ring = new RingAnalyzer(mol);

            double kM = 0;

            for (int idx = 0; idx < KmFragments.length; idx++) {

                double curVal = Descriptor.MISSING_VALUE;

                // Fragment 1
                if ((idx + 1) == 1)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 2
                if ((idx + 1) == 2)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 3
                if ((idx + 1) == 3)
                    curVal = CurFG.GetByName("nROH").getValue();

                // Fragment 4
                if ((idx + 1) == 4)
                    curVal = CurFG.GetByName("nArOH").getValue();

                // Fragment 5
                if ((idx + 1) == 5)
                    curVal = CurFG.GetByName("nRCOOH").getValue();

                // Fragment 6
                if ((idx + 1) == 6)
                    curVal = CurFG.GetByName("nRCHO").getValue() + CurFG.GetByName("nArCHO").getValue();

                // Fragment 7
                if ((idx + 1) == 7)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 8
                if ((idx + 1) == 8)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 9
                if ((idx + 1) == 9)
                    curVal = CurFG.GetByName("n135-Triazines").getValue();

                // Fragment 10
                if ((idx + 1) == 10)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 11
                if ((idx + 1) == 11)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 12
                if ((idx + 1) == 12)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 13
                if ((idx + 1) == 13)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 14
                if ((idx + 1) == 14)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 15
                if ((idx + 1) == 15)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 16
                if ((idx + 1) == 16)
                    curVal = CurFG.GetByName("nArNO2").getValue();

                // Fragment 17
                if ((idx + 1) == 17)
                    curVal = CurFG.GetByName("nRNH2").getValue() + CurFG.GetByName("nRNHR").getValue();

                // Fragment 18
                if ((idx + 1) == 18)
                    curVal = CurFG.GetByName("nArNH2").getValue() + CurFG.GetByName("nArNHR").getValue();

                // Fragment 19
                if ((idx + 1) == 19)
                    curVal = CurFG.GetByName("nRCN").getValue() + CurFG.GetByName("nArCN").getValue();

                // Fragment 20
                if ((idx + 1) == 20)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 21
                if ((idx + 1) == 21)
                    if (Queries[idx].matches(molStructure)) {
                        if ((Queries[idx].matchAll(molStructure).countUnique() == mol.GetStructure().getAtomCount())) {
                            Query_benzene.matches(molStructure);
                            if (Query_benzene.matchAll(molStructure).countUnique() >= 4)
                                curVal = 1;
                            else
                                curVal = 0;
                        } else
                            curVal = 0;
                    } else curVal = 0;

                // Fragment 22
                if ((idx + 1) == 22)
                    curVal = CurFG.GetByName("nPyridines").getValue();

                // Fragment 23
                if ((idx + 1) == 23)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 24
                if ((idx + 1) == 24)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 25
                if ((idx + 1) == 25)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 26
                if ((idx + 1) == 26)
                    curVal = CurFG.GetByName("nRNR2").getValue() + CurFG.GetByName("nArNR2").getValue();

                // Fragment 27
                if ((idx + 1) == 27)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 28
                if ((idx + 1) == 28)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 29
                if ((idx + 1) == 29)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 30
                if ((idx + 1) == 30)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 31
                if ((idx + 1) == 31) {
                    curVal = 0;
                    for (int k = 0; k < Queries_f_31.length; k++)
                        if (Queries_f_31[k].matches(molStructure))
                            curVal += Queries_f_31[k].matchAll(molStructure).countUnique();
                    if ((curVal > 0) && (curVal <= 3) && (ring.nRings <= 3))
                        curVal = 1;
                    else
                        curVal = 0;
                }

                // Fragment 32
                if ((idx + 1) == 32)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 33
                if ((idx + 1) == 33)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 34
                if ((idx + 1) == 34)
                    curVal = CurFG.GetByName("nCONN").getValue();

                // Fragment 35
                if ((idx + 1) == 35)
                    curVal = CurFG.GetByName("nFuranes").getValue() + CurFG.GetByName("nThiophenes").getValue();

                // Fragment 36
                if ((idx + 1) == 36)
                    curVal = CurFG.GetByName("nTriazoles").getValue();

                // Fragment 37
                if ((idx + 1) == 37)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 38
                if ((idx + 1) == 38)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 39
                if ((idx + 1) == 39)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 40
                if ((idx + 1) == 40)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 41
                if ((idx + 1) == 41)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 42
                if ((idx + 1) == 42)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 43
                if ((idx + 1) == 43)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 44
                if ((idx + 1) == 44)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 45
                if ((idx + 1) == 45)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 46
                if ((idx + 1) == 46)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 47
                if ((idx + 1) == 47)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique() * 3; // error in episuite?
                    else curVal = 0;

                // Fragment 48
                if ((idx + 1) == 48)
                    curVal = CurFG.GetByName("nThiazoles").getValue();

                // Fragment 49
                if ((idx + 1) == 49)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 50
                if ((idx + 1) == 50)
                    curVal = ring.nFusedAliphatic;

                // Fragment 51
                if ((idx + 1) == 51)
                    curVal = ring.nFusedBenzene;

                // Fragment 52
                if ((idx + 1) == 52)
                    if (ring.FourFusedAromatic) curVal = 1;
                    else curVal = 0;

                // Fragment 53
                if ((idx + 1) == 53)
                    if (ring.FourFusedAliphatic) curVal = 1;
                    else curVal = 0;

                // Fragment 54
                if ((idx + 1) == 54)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 55
                if ((idx + 1) == 55)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 56
                if ((idx + 1) == 56)
                    if (Queries[idx].matches(molStructure))
                        curVal = Queries[idx].matchAll(molStructure).countUnique();
                    else curVal = 0;

                // Fragment 57
                if ((idx + 1) == 57)
                    if (Queries[idx].matches(molStructure)) {
                        int nBiphenyls = Queries[idx].matchAll(molStructure).countUnique();
                        curVal = nBiphenyls;

                        // correction for benzenes (frag 54) previously found
                        double nBnz = GetByName("kM_f_54").getValue() - (nBiphenyls * 2);
                        nBnz = (nBnz > 0) ? nBnz : 0;
                        SetByName("kM_f_54", nBnz);
                        kM -= (nBiphenyls * 2) * (double) KmFragments[54 - 1][2];
                    } else curVal = 0;


                // finally set descriptor for fragment and update kM value
                SetByName("kM_f_" + (idx + 1), curVal);
                kM += curVal * (double) KmFragments[idx][2];
            }

            SetByName("kM", kM);

        } catch (Exception ex) {
            SetAllValues(Descriptor.MISSING_VALUE);
        }

    }


    /**
     * Clones the actual descriptor block
     *
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    @Override
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        KmFactor block = new KmFactor();
        block.CloneDetailsFrom(this);
        return block;
    }

}
