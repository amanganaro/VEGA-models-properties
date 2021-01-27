package insilico.guppy_knn.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.AtomCenteredFragments;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;

@Slf4j
public class ACF extends DescriptorBlock {

    private static final long serialVersionUID = 1L;


    private static final String BlockName = "Atom Centered Fragments";

    // Names and description of ACF

    private static String[] ACFNames = {
            "U-000", "C-001", "C-002", "C-003", "C-004", "C-005",
            "C-006", "C-007", "C-008", "C-009", "C-010", "C-011", "C-012",
            "C-013", "C-014", "C-015", "C-016", "C-017", "C-018", "C-019",
            "C-020", "C-021", "C-022", "C-023", "C-024", "C-025", "C-026",
            "C-027", "C-028", "C-029", "C-030", "C-031", "C-032", "C-033",
            "C-034", "C-035", "C-036", "C-037", "C-038", "C-039", "C-040",
            "C-041", "C-042", "C-043", "C-044", "U-045", "H-046", "H-047",
            "H-048", "H-049", "H-050", "H-051", "H-052", "H-053", "H-054",
            "H-055", "O-056", "O-057", "O-058", "O-059", "O-060", "O-061",
            "O-062", "O-063", "Se-064", "Se-065", "N-066", "N-067", "N-068",
            "N-069", "N-070", "N-071", "N-072", "N-073", "N-074", "N-075",
            "N-076", "N-077", "N-078", "N-079", "U-080", "F-081", "F-082",
            "F-083", "F-084", "F-085", "Cl-086", "Cl-087", "Cl-088", "Cl-089",
            "Cl-090", "Br-091", "Br-092", "Br-093", "Br-094", "Br-095",
            "I-096", "I-097", "I-098", "I-099", "I-100", "F-101", "Cl-102",
            "Br-103", "I-104", "U-105", "S-106", "S-107", "S-108", "S-109",
            "S-110", "Si-111",  "B-112", "U-113", "U-114", "P-115", "P-116",
            "P-117", "P-118", "P-119", "P-120" };

    private static String[] ACFDescription = {
            "undefined",
            "CH3R / CH4",
            "CH2R2",
            "CHR3",
            "CR4",
            "CH3X",
            "CH2RX",
            "CH2X2",
            "CHR2X",
            "CHRX2",
            "CHX3",
            "CR3X",
            "CR2X2",
            "CRX3",
            "CX4",
            " =CH2",
            " =CHR",
            " =CR2",
            " =CHX",
            " =CRX",
            " =CX2",
            "#CH",
            "#CR / R=C=R",
            "#CX",
            "R--CH--R",
            "R--CR--R",
            "R--CX--R",
            "R--CH--X",
            "R--CR--X",
            "R--CX--X",
            "X--CH--X",
            "X--CR--X",
            "X--CX--X",
            "R--CH..X",
            "R--CR..X",
            "R--CX..X",
            "Al-CH=X",
            "Ar-CH=X",
            "Al-C(=X)-Al",
            "Ar-C(=X)-R",
            "R-C(=X)-X / R-C#X / X=C=X",
            "X-C(=X)-X",
            "X--CH..X",
            "X--CR..X",
            "X--CX..X",
            "undefined",
            "H attached to C0(sp3) no X attached to next C",
            "H attached to C1(sp3)/C0(sp2)",
            "H attached to C2(sp3)/C1(sp2)/C0(sp)",
            "H attached to C3(sp3)/C2(sp2)/C3(sp2)/C3(sp)",
            "H attached to heteroatom",
            "H attached to alpha-C",
            "H attached to C0(sp3) with 1X attached to next C",
            "H attached to C0(sp3) with 2X attached to next C",
            "H attached to C0(sp3) with 3X attached to next C",
            "H attached to C0(sp3) with 4X attached to next C",
            "alcohol",
            "phenol / enol / carboxyl OH",
            "=O",
            "Al-O-Al",
            "Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X",
            "O--",
            "O- (negatively charged)",
            "R-O-O-R",
            "Any-Se-Any",
            "=Se",
            "Al-NH2",
            "Al2-NH",
            "Al3-N",
            "Ar-NH2 / X-NH2",
            "Ar-NH-Al",
            "Ar-NAl2",
            "RCO-N< / >N-X=X",
            "Ar2NH / Ar3N / Ar2N-Al / R..N..R",
            "R#N / R=N-",
            "R--N--R / R--N--X",
            "Ar-NO2 / R--N(--R)--O / RO-NO",
            "Al-NO2",
            "Ar-N=X / X-N=X",
            "N+ (positively charged)",
            "undefined",
            "F attached to C1(sp3)",
            "F attached to C2(sp3)",
            "F attached to C3(sp3)",
            "F attached to C1(sp2)",
            "F attached to C2(sp2)-C4(sp2)/C1(sp)/C4(sp3)/X",
            "Cl attached to C1(sp3)",
            "Cl attached to C2(sp3)",
            "Cl attached to C3(sp3)",
            "Cl attached to C1(sp2)",
            "Cl attached to C2(sp2)-C4(sp2)/C1(sp)/C4(sp3)/X",
            "Br attached to C1(sp3)",
            "Br attached to C2(sp3)",
            "Br attached to C3(sp3)",
            "Br attached to C1(sp2)",
            "Br attached to C2(sp2)-C4(sp2)/C1(sp)/C4(sp3)/X",
            "I attached to C1(sp3)",
            "I attached to C2(sp3)",
            "I attached to C3(sp3)",
            "I attached to C1(sp2)",
            "I attached to C2(sp2)-C4(sp2)/C1(sp)/C4(sp3)/X",
            "fluoride ion",
            "chloride ion",
            "bromide ion",
            "iodide ion",
            "undefined",
            "R-SH",
            "R2S / RS-SR",
            "R=S",
            "R-SO-R",
            "R-SO2-R",
            ">Si<",
            ">B- as in boranes",
            "undefined",
            "undefined",
            "P ylids",
            "R3-P=X",
            "X3-P=X (phosphate)",
            "PX3 (phosphite)",
            "PR3 (phosphine)",
            "C-P(X)2=X (phosphonate)"
    };

    // Private objects used inside the class only during calculation
    private IAtomContainer CurMol;
    private int nSK;
    private double[][] ConnAugMatrix;
    private boolean[] AtomAromatic;

    /**
     * Constructor.
     */
    public ACF() {
        super();
//        this.Name = AtomCenteredFragments.BlockName;
    }


    /**
     * Calculate descriptors for the given molecule.
     *
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {

        // Generate/clears descriptors
        GenerateDescriptors();

        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        // Gets matrices
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        int ACFSize = this.GetSize();

        int[] FragCount = new int[ACFSize];
        for (int i=0; i<ACFSize; i++)
            FragCount[i] = 0;

        // Inits basic data
        nSK = CurMol.getAtomCount();
        AtomAromatic = new boolean[nSK];
        for (int i=0; i<nSK; i++)
            AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

        // Searches fragments for each atom
        for (int i=0; i<nSK; i++) {
            CheckAtomFragments(i, FragCount);
        }

        // Fills descriptors DescList
        for (int i=0; i<ACFSize; i++)
            SetByIndex(i, FragCount[i]);

    }

    @Override
    protected void GenerateDescriptors() {
        DescList.clear();
        int ACFSize = ACFNames.length;
        for (int i=0; i<ACFSize; i++)
            Add(ACFNames[i], ACFDescription[i]);
        SetAllValues(Descriptor.MISSING_VALUE);
    }

    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    @Override
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        AtomCenteredFragments block = new AtomCenteredFragments();
//        block.CloneDetailsFrom(this);
        return block;
    }

    /**
     * Returns if an atomic number corresponds to an electronegative
     * atom (O, N, S, P, B, Si, Se or Halogens)
     *
     * @param AtomicNumber atomic number to be checked
     * @return true if it is electronegative
     */
    private boolean IsAtomElectronegative(int AtomicNumber) {
        boolean ans = false;
        if ((AtomicNumber==7)||(AtomicNumber==8)||(AtomicNumber==15)||
                (AtomicNumber==16)||(AtomicNumber==34)||(AtomicNumber==9)||
                (AtomicNumber==5)||(AtomicNumber==14)||
                (AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53)) {
            // O, N, S, P, B, Si, Se or Halogens (F, Cl, Br, I)
            ans = true;
        }
        return ans;
    }

    /**
     * Checks the atom number At from the CurMol (private object in this class)
     * and assigns it to an ACF, updating the counter array given as
     * parameter.
     *
     * @param At number of the atom to be checked
     * @param FragCount counter array to be updated
     */
    private void CheckAtomFragments(int At, int[] FragCount) {

        int VD=0, nH=0, Charge=0;
        int nSingle=0, nDouble=0, nTriple=0, nArom=0;

        int sX=0, dX=0, tX=0, aX=0, asX=0;
        int sR=0, dR=0, tR=0, aR=0;
        int sAr=0, dAr=0, aAr=0;
        int sAl=0, dAl=0, aAl=0;

        // note: asX is aromatic single bond like in pyrrole:
        // cn[H]c, csc, coc

        int C_OxiNumber=0, C_Hybridazion=0, C_CX=0;

        for (int j=0; j<nSK; j++) {
            if (j==At)
                continue;
            if (ConnAugMatrix[At][j]>0) {
                VD++;
                int Z = (int)ConnAugMatrix[j][j];
                double b = ConnAugMatrix[At][j];

                if (b==1) {
                    nSingle++;
                    if (IsAtomElectronegative(Z))
                        sX++;
                    if (Z==6)
                        sR++;
                    if (AtomAromatic[j])
                        sAr++;
                    else
                    if (Z==6) sAl++;
                }
                if (b==2) {
                    nDouble++;
                    if (IsAtomElectronegative(Z))
                        dX++;
                    if (Z==6)
                        dR++;
                    if (AtomAromatic[j])
                        dAr++;
                    else
                    if (Z==6) dAl++;
                }
                if (b==3) {
                    nTriple++;
                    if (IsAtomElectronegative(Z))
                        tX++;
                    if (Z==6)
                        tR++;
                }
                if (b==1.5) {
                    nArom++;
                    if (IsAtomElectronegative(Z)) {

                        // checks if is a pyrrole-like aromatic single bond
                        int elNegVD=0, elNegH=0, elNegCharge;
                        for (int k=0; k<nSK; k++)
                            if (ConnAugMatrix[j][k]>0) {
                                if (j==k) continue;
                                elNegVD++;
                            }
                        try {
                            elNegH = CurMol.getAtom(j).getImplicitHydrogenCount();
                        } catch (Exception e) {
                            elNegH = 0;
                        }
                        try {
                            elNegCharge = CurMol.getAtom(j).getFormalCharge();
                        } catch (Exception e) {
                            elNegCharge = 0;
                        }
                        elNegVD += elNegH - elNegCharge;

                        boolean IsPyrroleLikeArom=false;

                        if ((Z==7) && (elNegVD==3))
                            IsPyrroleLikeArom = true;
                        if (Z==8)
                            IsPyrroleLikeArom = true;
                        if ((Z==16) && (elNegVD==2) )
                            IsPyrroleLikeArom = true;

                        if (IsPyrroleLikeArom)
                            asX++;
                        else
                            aX++;

                    }
                    if (Z==6)
                        aR++;
                    if (AtomAromatic[j])
                        aAr++;
                    else
                    if (Z==6) aAl++;
                }
            }
        }

        // counts H
        try {
            nH = CurMol.getAtom(At).getImplicitHydrogenCount();
        } catch (Exception e) {
            nH = 0;
        }

        // formal charge
        try {
            Charge = CurMol.getAtom(At).getFormalCharge();
        } catch (Exception e) {
            Charge = 0;
        }

        // If Carbon, calculates oxidation number and hybridization
        if (ConnAugMatrix[At][At] == 6) {

            int c_VD=0;
            for (int j=0; j<nSK; j++) {
                if (j==At) continue;
                if (ConnAugMatrix[j][At]>0) {
                    c_VD++;
                    if (ConnAugMatrix[j][j] == 6) {
                        // search for -C-X
                        for (int k=0; k<nSK; k++) {
                            if (k==j) continue;
                            if ((ConnAugMatrix[k][j]>0) && (IsAtomElectronegative((int)ConnAugMatrix[k][k])))
                                C_CX++;
                        }
                    }
                }
            }

            C_OxiNumber += (sX) + (dX * 2) + (tX * 3);
            C_OxiNumber += (asX);
            if (aX > 1)
                C_OxiNumber += (aX * 1.5);
            else
                C_OxiNumber += (aX * 2);

            C_Hybridazion = c_VD + nH - 1;
            if (((int)C_OxiNumber)!=C_OxiNumber)
                C_OxiNumber = (int)C_OxiNumber + 1;
        }


        //// Search for proper fragment


        //// Hydrogen fragments ////////////////////////////////////////////////

        if (nH > 0) {

            if (ConnAugMatrix[At][At] == 6) {

                boolean IsAlphaCarbon = false;

                // Checks for alpha carbon

                if ((nSingle > 0) && (nDouble==0) && (nTriple==0) && (nArom==0)) {

                    for (int j=0; j<nSK; j++) {
                        if (j==At)
                            continue;

                        // -C
                        if (ConnAugMatrix[At][j] > 0) {

                            if ((ConnAugMatrix[At][j]==1) && (ConnAugMatrix[j][j]==6)) {

                                int nCdX=0, nCtX=0, nCaX=0;

                                for (int k=0; k<nSK; k++) {
                                    if (k==j)
                                        continue;
                                    int Z = (int) ConnAugMatrix[k][k];
                                    if ((ConnAugMatrix[j][k]>0) && (IsAtomElectronegative(Z))) {
                                        if (ConnAugMatrix[j][k]==2)
                                            nCdX++;
                                        if (ConnAugMatrix[j][k]==3)
                                            nCtX++;
                                        if (ConnAugMatrix[j][k]==1.5)
                                            nCaX++;
                                    }
                                }

                                if ((nCdX + nCtX + nCaX) == 1)
                                    IsAlphaCarbon = true;

                            } else {

                                IsAlphaCarbon = false; break;

                            }
                        }
                    }
                }

                if (IsAlphaCarbon) {

                    FragCount[51]+=nH;

                } else {

                    // C0sp3 (no X attached to next C)
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==0))
                        FragCount[46]+=nH;

                    // C1sp3, C0sp2
                    if ( ((C_OxiNumber==1) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==2)) )
                        FragCount[47]+=nH;

                    // C2sp3, C1sp2, C0sp
                    if ( ((C_OxiNumber==2) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==1) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==1)) )
                        FragCount[48]+=nH;

                    // C3sp3, C2sp2, C2sp2, C3sp
                    if ( ((C_OxiNumber==3) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==2) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==1)) )
                        FragCount[49]+=nH;

                    // C0sp3 with 1 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==1))
                        FragCount[52]+=nH;

                    // C0sp3 with 2 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==2))
                        FragCount[53]+=nH;

                    // C0sp3 with 3 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==3))
                        FragCount[54]+=nH;

                    // C0sp3 with 4 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==4))
                        FragCount[55]+=nH;

                }

            } else {

                // H to heteroatom
                FragCount[50]+=nH;

            }

        }



        //// Halogen fragments /////////////////////////////////////////////////

        // Search for halogens attached to current atom
        for (int j=0; j<nSK; j++) {
            if (j==At) continue;
            if (ConnAugMatrix[j][At]>0) {

                // F
                if (ConnAugMatrix[j][j] == 9) {

                    if (ConnAugMatrix[At][At] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragCount[81]++; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragCount[82]++; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragCount[83]++; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragCount[84]++; continue; }

                        // other cases fall into this fragment
                        FragCount[85]++; continue;

                    } else {

                        // attached to heteroatom
                        FragCount[85]++; continue;

                    }
                }

                // Cl
                if (ConnAugMatrix[j][j] == 17) {

                    if (ConnAugMatrix[At][At] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragCount[86]++; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragCount[87]++; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragCount[88]++; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragCount[89]++; continue; }

                        // other cases fall into this fragment
                        FragCount[90]++; continue;

                    } else {

                        // attached to heteroatom
                        FragCount[90]++; continue;

                    }
                }

                // Br
                if (ConnAugMatrix[j][j] == 35) {

                    if (ConnAugMatrix[At][At] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragCount[91]++; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragCount[92]++; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragCount[93]++; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragCount[94]++; continue; }

                        // other cases fall into this fragment
                        FragCount[95]++; continue;

                    } else {

                        // attached to heteroatom
                        FragCount[95]++; continue;

                    }
                }

                // I
                if (ConnAugMatrix[j][j] == 53) {

                    if (ConnAugMatrix[At][At] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragCount[96]++; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragCount[97]++; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragCount[98]++; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragCount[99]++; continue; }

                        // other cases fall into this fragment
                        FragCount[100]++;

                    } else {

                        // attached to heteroatom
                        FragCount[100]++;

                    }
                }

            }
        }



        //// Atom: C ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 6) {

            if ( ((nH==3) && (sR==1)) || (nH==4) )
            { FragCount[1]++; return; }

            if ((nH==2) && (sR==2))
            { FragCount[2]++; return; }

            if ((nH==1) && (sR==3))
            { FragCount[3]++; return; }

            if ((nH==0) && (sR==4))
            { FragCount[4]++; return; }

            if ((nH==3) && (sX==1))
            { FragCount[5]++; return; }

            if ((nH==2) && (sR==1) && (sX==1))
            { FragCount[6]++; return; }

            if ((nH==2) && (sX==2))
            { FragCount[7]++; return; }

            if ((nH==1) && (sR==2) && (sX==1))
            { FragCount[8]++; return; }

            if ((nH==1) && (sR==1) && (sX==2))
            { FragCount[9]++; return; }

            if ((nH==1) && (sX==3))
            { FragCount[10]++; return; }

            if ((nH==0) && (sR==3) && (sX==1))
            { FragCount[11]++; return; }

            if ((nH==0) && (sR==2) && (sX==2))
            { FragCount[12]++; return; }

            if ((nH==0) && (sR==1) && (sX==3))
            { FragCount[13]++; return; }

            if ((nH==0) && (sX==4))
            { FragCount[14]++; return; }

            if ((nH==2) && (dR==1))
            { FragCount[15]++; return; }

            if ((nH==1) && (dR==1) && (sR==1))
            { FragCount[16]++; return; }

            if ((nH==0) && (dR==1) && (sR==2))
            { FragCount[17]++; return; }

            if ((nH==1) && (dR==1) && (sX==1))
            { FragCount[18]++; return; }

            if ((nH==0) && (dR==1) && (sX==1) && (sR==1))
            { FragCount[19]++; return; }

            if ((nH==0) && (dR==1) && (sX==2))
            { FragCount[20]++; return; }

            if ((nH==1) && (tR==1))
            { FragCount[21]++; return; }

            if ( ((nH==0) && (tR==1) && (sR==1)) ||
                    ((nH==0) && (dR==2)) )
            { FragCount[22]++; return; }

            if ((nH==0) && (tR==1) && (sX==1))
            { FragCount[23]++; return; }

            if ((nH==1) && (VD==2) && (aR==2))
            { FragCount[24]++; return; }

            if ((nH==0) && (VD==3) && (aR>=2) && ((sR==1)||(aR==3)))
            { FragCount[25]++; return; }

            if ((nH==0) && (VD==3) && (aR==2) && (sX==1))
            { FragCount[26]++; return; }

            if ((nH==1) && (VD==2) && (aR==1) && (aX==1))
            { FragCount[27]++; return; }

            if ((nH==0) && (VD==3) && (aR>=1) && (aX==1) && ((sR==1)||(aR==2)))
            { FragCount[28]++; return; }

            if ((nH==0) && (VD==3) && (aR==1) && (aX==1) && (sX==1))
            { FragCount[29]++; return; }

            if ((nH==1) && (VD==2) && (aX==2))
            { FragCount[30]++; return; }

            if ((nH==0) && (VD==3) && (aX==2) && ((sR==1)||(aR==1)))
            { FragCount[31]++; return; }

            if ((nH==0) && (VD==3) && (aX==2) && (sX==1))
            { FragCount[32]++; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aR==1) && (asX==1))
            { FragCount[33]++; return; }

            if ((nH==0) && (VD==3) && (((sR==1) && (aR==1))||(aR==2)) && (asX==1))
            { FragCount[34]++; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aR==1) && (asX==1))
            { FragCount[35]++; return; }

            ///////////////////////////////////////////////////////


            if ((nH==1) && (dX==1) && (sAl==1))
            { FragCount[36]++; return; }

            if ((nH==1) && (dX==1) && (sAr==1))
            { FragCount[37]++; return; }

            if ((nH==0) && (dX==1) && (sAl==2))
            { FragCount[38]++; return; }

            if ((nH==0) && (dX==1) && (sAr==1) && (sR>=1) &&(sX==0))
            { FragCount[39]++; return; }

            if (((nH==0) && (dX==1) && (sR==1) && (sX==1)) ||
                    ((nH==0) && (tX==1) && (sR==1)) ||
                    ((nH==0) && (dX==2)) )
            { FragCount[40]++; return; }

            if ((nH==0) && (dX==1) && (sX==2))
            { FragCount[41]++; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aX==1) && (asX==1))
            { FragCount[42]++; return; }

            if ((nH==0) && (VD==3) && ((sR==1)||(aR==1)) && (aX==1) && (asX==1))
            { FragCount[43]++; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aX==1) && (asX==1))
            { FragCount[44]++; return; }

            ///////////////////////////////////////////////////////


            return;
        }



        //// Atom: O ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 8) {

            // Checks for particular O fragments linked to C or N

            for (int j=0; j<nSK; j++) {
                if (j==At)
                    continue;

                // Bound to C

                if ((ConnAugMatrix[At][j]>0) && (ConnAugMatrix[j][j]==6)) {
                    boolean c_Arom = AtomAromatic[j];
                    int c_VD = 0, c_dR=0, c_dO=0, c_dX=0;
                    for (int k=0; k<nSK; k++) {
                        if (k==j)
                            continue;
                        if (ConnAugMatrix[j][k]>0)  {
                            c_VD++;
                            if (ConnAugMatrix[j][k] == 2) {
                                if (IsAtomElectronegative((int)ConnAugMatrix[k][k]))
                                    c_dX++;
                                if (ConnAugMatrix[k][k] == 8)
                                    c_dO++;
                                else
                                    c_dR++;
                            }
                        }
                    }


                    // R-O-C=X (for fragment O-060

                    if (((VD+nH)==2) && (nSingle==2) && (c_dX>0))
                    { FragCount[60]++; return; }


                    // Particular OH groups

                    if (nH == 1) {
                        // Enol
                        if ((c_VD==3) && (!c_Arom) && (c_dR==1))
                        { FragCount[57]++; return; }

                        // Phenol
                        if (c_Arom)
                        { FragCount[57]++; return; }

                        // Carboxyl
                        if ((!c_Arom) && (c_VD==3) && (c_dO==1))
                        { FragCount[57]++; return; }

                    }

                }


                // Bound to N

                if ((ConnAugMatrix[At][j]>0) && (ConnAugMatrix[j][j]==7)) {
                    int n_VD = 0, n_dO=0, n_sOminus=0;
                    for (int k=0; k<nSK; k++) {
                        if (k==j)
                            continue;
                        if (ConnAugMatrix[j][k]>0)  {
                            n_VD++;
                            if (ConnAugMatrix[j][k] == 2) {
                                if (ConnAugMatrix[k][k] == 8)
                                    n_dO++;
                            }
                            if (ConnAugMatrix[j][k] == 1) {
                                if (ConnAugMatrix[k][k] == 8) {
                                    int o_VD=0;
                                    for (int z=0; z<nSK; z++) {
                                        if (z==k) continue;
                                        if (ConnAugMatrix[k][z]>0) o_VD++;
                                    }
                                    if (o_VD == 1) {
                                        int nH_O = 0;
                                        try {
                                            nH_O = CurMol.getAtom(k).getImplicitHydrogenCount();
                                        } catch (Exception e) {
                                            nH_O = 0;
                                        }
                                        if (nH == 0)
                                            n_sOminus++;
                                    }
                                }
                            }
                        }
                    }


                    // Nitro group
                    // in NO2 both oxygens are seen as O--

                    if ((VD==1) && (nH==0) && (nDouble==1) && (n_VD==3) && (n_sOminus==1))
                    { FragCount[61]++; return; }

                    if ((VD==1) && (nH==0) && (n_VD==3) && (n_dO==1))
                    { FragCount[61]++; return; }

                }

            }

            // generic alcohol (even if not bound to C)
            if ((nSingle==1) && (nH==1))
            { FragCount[56]++; return; }


            if ((VD==1) && (nDouble==1))
            { FragCount[58]++; return; }

            if ((VD==2) && (nSingle==2)) {
                for (int j=0; j<nSK; j++) {
                    if (j==At)
                        continue;
                    if ((ConnAugMatrix[At][j]==1) && (ConnAugMatrix[j][j]==8)) {
                        // Found an oxygen, R-O-O
                        FragCount[63]++;
                        return;
                    }
                }
            }

            if (((VD==2) && (sAl==1) && (sAr==1)) ||
                    ((VD==2) && (sAr==2)) ||
                    ((VD==2) && (aR==2)))
            { FragCount[60]++; return; }

            if ((VD==2) && (nSingle==2) && (sAr==0))
            { FragCount[59]++; return; }

            if (Charge == -2)
            { FragCount[61]++; return; }

            if (Charge == -1)
            { FragCount[62]++; return; }

            return;
        }



        //// Atom: Se ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 34) {

            if ((VD==2) && (nSingle==2))
            { FragCount[64]++; return; }

            if ((VD==1) && (nDouble==1))
            { FragCount[65]++; return; }

            return;
        }



        //// Atom: N ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 7) {

            // Checks for particular aromatic form
            boolean AromPyridineLike=false, AromPyrroleLike=false;
            if (nArom >= 2) {
                if ((VD + nH - Charge) == 2 )
                    AromPyridineLike = true;
                if ((VD + nH - Charge) == 3 )
                    AromPyrroleLike = true;
            }

            // Checks for particular N fragments linked to C, O or X

            int nO=0, sOCR=0, sXdX=0;

            for (int j=0; j<nSK; j++) {
                if (j==At)
                    continue;
                if (ConnAugMatrix[At][j]>0) {
                    int Z = (int)ConnAugMatrix[j][j];

                    // O
                    if (Z==8) {
                        nO++;
                        for (int k=0; k<nSK; k++) {
                            if ((k==j) || (k==At)) continue;
                            if ((ConnAugMatrix[k][j]>0) && (ConnAugMatrix[k][k]==6)) {
                                int c_VD=0;
                                for (int z=0; z<nSK; z++) {
                                    if (z==k) continue;
                                    if (ConnAugMatrix[z][k]>0) c_VD++;
                                }
                                if (c_VD == 2)
                                    sOCR++;
                            }
                        }
                    }

                    // Electronegative
                    // NOTE: to be consistent with Dragon, also R is used
                    // instead of only electronegative atoms
                    if ((IsAtomElectronegative(Z)) || (Z==6) ) {
                        int x_VD=0, x_dX=0;
                        for (int k=0; k<nSK; k++) {
                            if ((k==j) || (k==At)) continue;
                            if (ConnAugMatrix[k][j]>0) {
                                x_VD++;
                                if (IsAtomElectronegative((int)ConnAugMatrix[k][k]))
                                    if (ConnAugMatrix[k][j] == 2)
                                        x_dX++;
                            }
                        }
//                        if ((x_VD==1) && (x_dX==1))
                        if (x_dX==1)
                            sXdX++;
                    }

                }
            }

            // fragments with possible [N+](=O)[O-] or R-O-N=O
            if (((VD==3) && (nH==0) && (sAr==1) && (nO==2)) ||
                    ((VD==3) && (nH==0) && (aR==2) && (aX==1)) ||
                    ((VD==2) && (nH==0) && (nO==2)))
//                || ((VD==2) && (nH==0) && (nO==2)))
            { FragCount[76]++; return; }

            if ((VD==3) && (nH==0) && (sAl==1) && (nO==2))
            { FragCount[77]++; return; }

            // N Charged +1
            if ((!AromPyrroleLike) && (!AromPyridineLike))
                if (Charge==1)
                { FragCount[79]++; return; }

            // fragment with particular groups
            if (!AromPyrroleLike)
                if (((VD==3) && (sOCR==1)) ||   // >N-OCR
                        (((VD+nH)==3) && (sXdX>0)))     // >N-X=X  NOTE: also >N-R=X for Dragon compatibility
                { FragCount[72]++; return; }


            if ((VD==1) && (nH==2) && (sAl==1))
            { FragCount[66]++; return; }

            if ((VD==2) && (nH==1) && (sAl==2))
            { FragCount[67]++; return; }

            if ((VD==3) && (nH==0) && (sAl==3))
            { FragCount[68]++; return; }

            if ((VD==1) && (nH==2) && ( (sAr==1) || (sX==1) ))
            { FragCount[69]++; return; }

            if ((VD==2) && (nH==1) && (sAl==2))
            { FragCount[70]++; return; }

            if ((VD==3) && (nH==0) && (sAl==2) && (sAr==1))
            { FragCount[71]++; return; }

            if (((VD==2) && (nH==1) && (sAr==2)) ||
                    ((VD==3) && (nH==0) && (sAr==3)) ||
                    ((VD==3) && (nH==0) && (sAr==2) && (sAl==1)) ||
                    ((AromPyrroleLike)))
            { FragCount[73]++; return; }

            if (((VD==1) && (nH==0) && (tR==1)) ||
                    ((VD==2) && (nH==0) && (dR==1) && (sR==1)))
            { FragCount[74]++; return; }

            if (((AromPyridineLike)) ) // ||
//                ((VD==2) && (aR==1) && (aX==1)))
            { FragCount[75]++; return; }

            if (((VD==2) && (nH==0) && (sAr==1) && (dX==1)) ||
                    ((VD==2) && (nH==0) && (sX==1) && (dX==1)))
            { FragCount[78]++; return; }

            return;
        }



        //// Halogen ions //////////////////////////////////////////////////////

        // F
        if ((ConnAugMatrix[At][At] == 9) && (VD==0))
        { FragCount[101]++; return; }

        // Cl
        if ((ConnAugMatrix[At][At] == 17) && (VD==0))
        { FragCount[102]++; return; }

        // Br
        if ((ConnAugMatrix[At][At] == 35) && (VD==0))
        { FragCount[103]++; return; }

        // I
        if ((ConnAugMatrix[At][At] == 53) && (VD==0))
        { FragCount[104]++; return; }



        //// Atom: S ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 16) {

            int nO=0;

            for (int j=0; j<nSK; j++) {
                if (j==At) continue;
                if (ConnAugMatrix[At][j]>0) {
                    if (ConnAugMatrix[j][j]==8)
                        nO++;
                }
            }

            if ((nH==1) && (sR==1))
            { FragCount[106]++; return; }

            if ((VD==2) && (nH==0))    // RSR is valid both for single and aromatic bonds (?)
            { FragCount[107]++; return; }

            if ((VD==1) && (nDouble==1))
            { FragCount[108]++; return; }

            if ((VD==3) && (nSingle==2) && (nO==1))
            { FragCount[109]++; return; }

            if ((VD==4) && (nSingle==2) && (nO>=2))
            { FragCount[110]++; return; }

            return;
        }



        //// Atom: Si //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 14) {

            if ((VD==4) && (nSingle==4))
            { FragCount[111]++; return; }

            return;

        }



        //// Atom: B //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 5) {

            if ((VD==3) && (nSingle==3))
            { FragCount[112]++; return; }

            return;

        }



        //// Atom: P //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 15) {

            if ((VD==4) && (Charge==1))
            { FragCount[115]++; return; }

            if ((VD==4) && (dX==1) && (sR==3))
            { FragCount[116]++; return; }

            if ((VD==4) && (dX==1) && (sX==3))
            { FragCount[117]++; return; }


            if ((VD==3) && (sX==3))
            { FragCount[118]++; return; }

            if ((VD==3) && (sR==3))
            { FragCount[119]++; return; }

            if ((VD==4) && (dX==1) && (sR==1) && (sX==2))
            { FragCount[120]++; return; }

        }
    }
    
}
