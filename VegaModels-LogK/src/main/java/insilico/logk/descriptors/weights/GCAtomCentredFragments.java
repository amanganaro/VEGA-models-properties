package insilico.logk.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;

/**
 * Implementation of the Ghose-Crippen ACF, allowing to assign each atom to the specific fragment id.
 *
 */
public class GCAtomCentredFragments {

    //
    // Name and description of each ACF, as:
    // Array index, Fragment Symbol, Fragment Description
    //
    public final static String[][] ACF_NAMES = {
            { "0", "U-000", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "1", "C-001", "CH3R / CH4"},
            { "2", "C-002", "CH2R2"},
            { "3", "C-003", "CHR3"},
            { "4", "C-004", "CR4"},
            { "5", "C-005", "CH3X"},
            { "6", "C-006", "CH2RX"},
            { "7", "C-007", "CH2X2"},
            { "8", "C-008", "CHR2X"},
            { "9", "C-009", "CHRX2"},
            { "10", "C-010", "CHX3"},
            { "11", "C-011", "CR3X"},
            { "12", "C-012", "CR2X2"},
            { "13", "C-013", "CRX3"},
            { "14", "C-014", "CX4"},
            { "15", "C-015", "=CH2"},
            { "16", "C-016", "=CHR"},
            { "17", "C-017", "=CR2"},
            { "18", "C-018", "=CHX"},
            { "19", "C-019", "=CRX"},
            { "20", "C-020", "=CX2"},
            { "21", "C-021", "#CH"},
            { "22", "C-022", "#CR / R=C=R"},
            { "23", "C-023", "#CX"},
            { "24", "C-024", "R--CH--R"},
            { "25", "C-025", "R--CR--R"},
            { "26", "C-026", "R--CX--R"},
            { "27", "C-027", "R--CH--X"},
            { "28", "C-028", "R--CR--X"},
            { "29", "C-029", "R--CX--X"},
            { "30", "C-030", "X--CH--X"},
            { "31", "C-031", "X--CR--X"},
            { "32", "C-032", "X--CX--X"},
            { "33", "C-033", "R--CH..X"},
            { "34", "C-034", "R--CR..X"},
            { "35", "C-035", "R--CX..X"},
            { "36", "C-036", "Al-CH=X"},
            { "37", "C-037", "Ar-CH=X"},
            { "38", "C-038", "Al-C(=X)-Al"},
            { "39", "C-039", "Ar-C(=X)-R"},
            { "40", "C-040", "R-C(=X)-X / R-C#X / X=C=X"},
            { "41", "C-041", "X-C(=X)-X"},
            { "42", "C-042", "X--CH..X"},
            { "43", "C-043", "X--CR..X"},
            { "44", "C-044", "X--CX..X"},
            { "45", "U-045", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "46", "H-046", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_46")},
            { "47", "H-047", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_47")},
            { "48", "H-048", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_48")},
            { "49", "H-049", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_49")},
            { "50", "H-050", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_50")},
            { "51", "H-051", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_51")},
            { "52", "H-052", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_52")},
            { "53", "H-053", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_53")},
            { "54", "H-054", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_54")},
            { "55", "H-055", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_55")},
            { "56", "O-056", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_56")},
            { "57", "O-057", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_57")},
            { "58", "O-058", "=O"},
            { "59", "O-059", "Al-O-Al"},
            { "60", "O-060", "Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X"},
            { "61", "O-061", "O--"},
            { "62", "O-062", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_62")},
            { "63", "O-063", "R-O-O-R"},
            { "64", "Se-064", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_64")},
            { "65", "Se-065", "=Se"},
            { "66", "N-066", "Al-NH2"},
            { "67", "N-067", "Al2-NH"},
            { "68", "N-068", "Al3-N"},
            { "69", "N-069", "Ar-NH2 / X-NH2"},
            { "70", "N-070", "Ar-NH-Al"},
            { "71", "N-071", "Ar-NAl2"},
            { "72", "N-072", "RCO-N< / >N-X=X"},
            { "73", "N-073", "Ar2NH / Ar3N / Ar2N-Al / R..N..R"},
            { "74", "N-074", "R#N / R=N-"},
            { "75", "N-075", "R--N--R / R--N--X"},
            { "76", "N-076", "Ar-NO2 / R--N(--R)--O / RO-NO"},
            { "77", "N-077", "Al-NO2"},
            { "78", "N-078", "Ar-N=X / X-N=X"},
            { "79", "N-079", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_79")},
            { "80", "U-080", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "81", "F-081", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_81")},
            { "82", "F-082", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_82")},
            { "83", "F-083", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_83")},
            { "84", "F-084", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_84")},
            { "85", "F-085", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_85")},
            { "86", "Cl-086", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_86")},
            { "87", "Cl-087", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_87")},
            { "88", "Cl-088", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_88")},
            { "89", "Cl-089", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_89")},
            { "90", "Cl-090", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_90")},
            { "91", "Br-091", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_91")},
            { "92", "Br-092", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_92")},
            { "93", "Br-093", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_93")},
            { "94", "Br-094", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_94")},
            { "95", "Br-095", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_95")},
            { "96", "I-096", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_96")},
            { "97", "I-097", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_97")},
            { "98", "I-098", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_98")},
            { "99", "I-099", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_99")},
            { "100", "I-100", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_100")},
            { "101", "F-101", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_101")},
            { "102", "Cl-102", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_102")},
            { "103", "Br-103", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_103")},
            { "104", "I-104", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_104")},
            { "105", "U-105", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "106", "S-106", "R-SH"},
            { "107", "S-107", "R2S / RS-SR"},
            { "108", "S-108", "R=S"},
            { "109", "S-109", "R-SO-R"},
            { "110", "S-110", "R-SO2-R"},
            { "111", "Si-111", ">Si<"},
            { "112", "B-112", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_112")},
            { "113", "U-113", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "114", "U-114", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_undefined")},
            { "115", "P-115", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_115")},
            { "116", "P-116", "R3-P=X"},
            { "117", "P-117", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_117")},
            { "118", "P-118", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_118")},
            { "119", "P-119", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_119")},
            { "120", "P-120", StringSelectorCore.getString("molecule_acf_ghosecrippen_name_120")}
    };


    //
    // Names of the undefined fragments in the ACF_NAMES list
    //
    public final static String[] ACF_UNDEFINED_NAMES = {
            "U-000", "U-045", "U-080", "U-105", "U-113", "U-114"
    };


    private final IAtomContainer CurMol;
    private final int nSK;
    private final double[][] ConnAugMatrix;
    private final boolean[] AtomAromatic;
    private final boolean ExplicitHydrogen;
    private final HashMap<Integer, Integer> NotMappedFragCount;

    private final int[] FragAtomId;


    public GCAtomCentredFragments(IAtomContainer Mol, boolean HasExplicitHydrogen) {

        // Init all variables
        this.CurMol = Mol;
        this.ExplicitHydrogen = HasExplicitHydrogen;
        nSK = Mol.getAtomCount();
        NotMappedFragCount = new HashMap<>();
        ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
        AtomAromatic = new boolean[nSK];
        for (int i=0; i<nSK; i++)
            AtomAromatic[i] = CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC);

        // Init array of ACF id, all to MV
        FragAtomId = new int[nSK];
        for (int i : FragAtomId)
            i = Descriptor.MISSING_VALUE;

        // Check all atoms and try to assign to an ACF
        for (int AtomIdx=0; AtomIdx<nSK; AtomIdx++)
            ProcessAtom(AtomIdx);

    }


    public int[] GetACF() {
        return this.FragAtomId;
    }

    /**
     * Provide an hashmap with all the fragments (id, count) that are not directly mapped in the molecule. This is
     * needed to retrieve H fragments when a H-depleted molecule is processed.
     *
     * @return
     */
    public HashMap<Integer, Integer> getNotMappedFragCount() { return this.NotMappedFragCount; }


    private void ProcessAtom(int AtomIndex) {

        IAtom CurAt = CurMol.getAtom(AtomIndex);

        // if H skip - it will be calculated from the attached atom
        if (ExplicitHydrogen)
            if (ConnAugMatrix[AtomIndex][AtomIndex] == 1)
                return;

        // if halogen skip - it will be calculated from the attached atom
        if (IsAtomHalogen(CurAt.getAtomicNumber()))
            return;


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
            if (j==AtomIndex)
                continue;
            if ( (ConnAugMatrix[AtomIndex][j]>0) && (ConnAugMatrix[j][j] != 1) ) {
                VD++;
                int Z = (int)ConnAugMatrix[j][j];
                double b = ConnAugMatrix[AtomIndex][j];

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
                            if ((ConnAugMatrix[j][k]>0) && (ConnAugMatrix[k][k]!=1)) {
                                if (j==k) continue;
                                elNegVD++;
                            }
                        try {
                            if (ExplicitHydrogen) {
                                for (IAtom connAt : CurMol.getConnectedAtomsList(CurMol.getAtom(j)))
                                    if (connAt.getAtomicNumber() == 1)
                                        elNegH++;
                            } else {
                                elNegH = CurMol.getAtom(j).getImplicitHydrogenCount();
                            }
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
                        if ((Z==8) && (elNegVD==2))
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
        if (ExplicitHydrogen) {
            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex) continue;
                if (ConnAugMatrix[j][AtomIndex] == 1)
                    if (ConnAugMatrix[j][j] == 1)
                        nH++;
            }
        } else {
            try {
                nH = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                nH = 0;
            }
        }

        // formal charge
        try {
            Charge = CurAt.getFormalCharge();
        } catch (Exception e) {
            Charge = 0;
        }

        // If Carbon, calculates oxidation number and hybridization
        if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

            int c_VD=0;
            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex) continue;
                if ( (ConnAugMatrix[j][AtomIndex]>0) && (ConnAugMatrix[j][j] != 1) ) {
                    c_VD++;
//                    if (IsAtomElectronegative((int)ConnMatrix[j][j])) {
//                        // Electronegative heteroatom, bond value added to oxidation
//                        double OxiAddVal = ConnMatrix[j][At];
//                        if (OxiAddVal==1.5) {
//                            // Checks for aromatic heteroatom with 2 pi electron
//                            int he_VD=0;
//                            for (int k=0; k<nSK; k++) {
//                                if (k==j) continue;
//                                if (ConnMatrix[j][k]>0) he_VD++;
//                            }
//                            int he_nH=0;
//                            try {
//                                he_nH = CurMol.getAtom(j).getHydrogenCount();
//                            } catch (Exception e) {  }
//                            he_VD += he_nH;
//                            if (he_VD == 3)
//                                OxiAddVal = 1;
//                        }
//                        C_OxiNumber += OxiAddVal;
//                    }
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

            int H_type = 0;

            if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

                boolean IsAlphaCarbon = false;

                // Checks for alpha carbon

                if ((nSingle > 0) && (nDouble==0) && (nTriple==0) && (nArom==0)) {

                    for (int j=0; j<nSK; j++) {
                        if (j==AtomIndex)
                            continue;

                        // -C
                        if (ConnAugMatrix[AtomIndex][j] > 0) {

                            if (ConnAugMatrix[j][j]==1)
                                continue;

                            if ((ConnAugMatrix[AtomIndex][j]==1) && (ConnAugMatrix[j][j]==6)) {

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

                    H_type=51;

                } else {

                    // C0sp3 (no X attached to next C)
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==0))
                        H_type=46;

                    // C1sp3, C0sp2
                    if ( ((C_OxiNumber==1) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==2)) )
                        H_type=47;

                    // C2sp3, C1sp2, C0sp
                    if ( ((C_OxiNumber==2) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==1) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==1)) )
                        H_type=48;

                    // C3sp3, C2sp2, C2sp2, C3sp
                    if ( ((C_OxiNumber==3) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==2) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==1)) )
                        H_type=49;

                    // C0sp3 with 1 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==1))
                        H_type=52;

                    // C0sp3 with 2 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==2))
                        H_type=53;

                    // C0sp3 with 3 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==3))
                        H_type=54;

                    // C0sp3 with 4 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==4))
                        H_type=55;

                }

            } else {

                // H to heteroatom
                H_type=50;

            }

            if (ExplicitHydrogen) {

                // Sets the id for H atoms as they are explicit
                for (int idxH=0; idxH<nSK; idxH++) {
                    if (idxH == AtomIndex) continue;
                    if ((ConnAugMatrix[AtomIndex][idxH] == 1) && (ConnAugMatrix[idxH][idxH] == 1) )
                        FragAtomId[idxH] = H_type;
                }

            } else {

                // Implicit H: counts the number of H and store the value in an hashmap, as these atoms
                // can not be mapped directly on the molecule
                int CountValue = nH;
                if (NotMappedFragCount.containsKey(H_type))
                    CountValue += NotMappedFragCount.get(H_type);
                NotMappedFragCount.put(H_type, CountValue);
            }
        }



        //// Halogen fragments /////////////////////////////////////////////////

        // Search for halogens attached to current atom
        for (int j=0; j<nSK; j++) {
            if (j==AtomIndex) continue;
            if (ConnAugMatrix[j][AtomIndex]>0) {

                // F
                if (ConnAugMatrix[j][j] == 9) {

                    if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragAtomId[j] = 81; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragAtomId[j] = 82; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragAtomId[j] = 83; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragAtomId[j] = 84; continue; }

                        // other cases fall into this fragment
                        FragAtomId[j] = 85; continue;

                    } else {

                        // attached to heteroatom
                        FragAtomId[j] = 85; continue;

                    }
                }

                // Cl
                if (ConnAugMatrix[j][j] == 17) {

                    if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragAtomId[j] = 86; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragAtomId[j] = 87; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragAtomId[j] = 88; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragAtomId[j] = 89; continue; }

                        // other cases fall into this fragment
                        FragAtomId[j] = 90; continue;

                    } else {

                        // attached to heteroatom
                        FragAtomId[j] = 90; continue;

                    }
                }

                // Br
                if (ConnAugMatrix[j][j] == 35) {

                    if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragAtomId[j] = 91; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragAtomId[j] = 92; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragAtomId[j] = 93; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragAtomId[j] = 94; continue; }

                        // other cases fall into this fragment
                        FragAtomId[j] = 95; continue;

                    } else {

                        // attached to heteroatom
                        FragAtomId[j] = 95; continue;

                    }
                }

                // I
                if (ConnAugMatrix[j][j] == 53) {

                    if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {
                        // attached to C1sp3
                        if ((C_OxiNumber==1) && (C_Hybridazion==3))
                        { FragAtomId[j] = 96; continue; }

                        // attached to C2sp3
                        if ((C_OxiNumber==2) && (C_Hybridazion==3))
                        { FragAtomId[j] = 97; continue; }

                        // attached to C3sp3
                        if ((C_OxiNumber==3) && (C_Hybridazion==3))
                        { FragAtomId[j] = 98; continue; }

                        // attached to C1sp2
                        if ((C_OxiNumber==1) && (C_Hybridazion==2))
                        { FragAtomId[j] = 99; continue; }

                        // other cases fall into this fragment
                        FragAtomId[j] = 100; continue;

                    } else {

                        // attached to heteroatom
                        FragAtomId[j] = 100; continue;

                    }
                }

            }
        }



        //// Atom: C ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 6) {

            if ( ((nH==3) && (sR==1)) || (nH==4) )
            { FragAtomId[AtomIndex] = 1; return; }

            if ((nH==2) && (sR==2))
            { FragAtomId[AtomIndex] = 2; return; }

            if ((nH==1) && (sR==3))
            { FragAtomId[AtomIndex] = 3; return; }

            if ((nH==0) && (sR==4))
            { FragAtomId[AtomIndex] = 4; return; }

            if ((nH==3) && (sX==1))
            { FragAtomId[AtomIndex] = 5; return; }

            if ((nH==2) && (sR==1) && (sX==1))
            {FragAtomId[AtomIndex] = 6; return; }

            if ((nH==2) && (sX==2))
            { FragAtomId[AtomIndex] = 7; return; }

            if ((nH==1) && (sR==2) && (sX==1))
            { FragAtomId[AtomIndex] = 8; return; }

            if ((nH==1) && (sR==1) && (sX==2))
            { FragAtomId[AtomIndex] = 9; return; }

            if ((nH==1) && (sX==3))
            { FragAtomId[AtomIndex] = 10; return; }

            if ((nH==0) && (sR==3) && (sX==1))
            { FragAtomId[AtomIndex] = 11; return; }

            if ((nH==0) && (sR==2) && (sX==2))
            { FragAtomId[AtomIndex] = 12; return; }

            if ((nH==0) && (sR==1) && (sX==3))
            { FragAtomId[AtomIndex] = 13; return; }

            if ((nH==0) && (sX==4))
            { FragAtomId[AtomIndex] = 14; return; }

            if ((nH==2) && (nDouble==1) && (dX==0))
            { FragAtomId[AtomIndex] = 15; return; }

            if ((nH==1) && (nDouble==1) && (sR==1) && (dX==0))
            { FragAtomId[AtomIndex] = 16; return; }

            if ((nH==0) && (nDouble==1) && (sR==2) && (dX==0))
            { FragAtomId[AtomIndex] = 17; return; }

            if ((nH==1) && (nDouble==1) && (sX==1) && (dX==0))
            { FragAtomId[AtomIndex] = 18; return; }

            if ((nH==0) && (nDouble==1) && (sX==1) && (sR==1) && (dX==0))
            { FragAtomId[AtomIndex] = 19; return; }

            if ((nH==0) && (nDouble==1) && (sX==2) && (dX==0))
            { FragAtomId[AtomIndex] = 20; return; }

            if ((nH==1) && (tR==1))
            { FragAtomId[AtomIndex] = 21; return; }

            if ( ((nH==0) && (tR==1) && (sR==1)) ||
                    ((nH==0) && (dR==2)) )
            { FragAtomId[AtomIndex] = 22; return; }

            if ((nH==0) && (tR==1) && (sX==1))
            { FragAtomId[AtomIndex] = 23; return; }

            if ((nH==1) && (VD==2) && (aR==2))
            { FragAtomId[AtomIndex] = 24; return; }

            if ((nH==0) && (VD==3) && (aR>=2) && ((sR==1)||(aR==3)))
            { FragAtomId[AtomIndex] = 25; return; }

            if ((nH==0) && (VD==3) && (aR==2) && (sX==1))
            { FragAtomId[AtomIndex] = 26; return; }

            if ((nH==1) && (VD==2) && (aR==1) && (aX==1))
            { FragAtomId[AtomIndex] = 27; return; }

            if ((nH==0) && (VD==3) && (aR>=1) && (aX==1) && ((sR==1)||(aR==2)))
            { FragAtomId[AtomIndex] = 28; return; }

            if ((nH==0) && (VD==3) && (aR==1) && (aX==1) && (sX==1))
            { FragAtomId[AtomIndex] = 29; return; }

            if ((nH==1) && (VD==2) && (aX==2))
            { FragAtomId[AtomIndex] = 30; return; }

            if ((nH==0) && (VD==3) && (aX==2) && ((sR==1)||(aR==1)))
            { FragAtomId[AtomIndex] = 31; return; }

            if ((nH==0) && (VD==3) && (aX==2) && (sX==1))
            { FragAtomId[AtomIndex] = 32; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aR==1) && (asX==1))
            { FragAtomId[AtomIndex] = 33; return; }

            if ((nH==0) && (VD==3) && (((sR==1) && (aR==1))||(aR==2)) && (asX==1))
            { FragAtomId[AtomIndex] = 34; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aR==1) && (asX==1))
            { FragAtomId[AtomIndex] = 35; return; }

            ///////////////////////////////////////////////////////


            if ((nH==1) && (dX==1) && (sAl==1))
            { FragAtomId[AtomIndex] = 36; return; }

            if ((nH==1) && (dX==1) && (sAr==1))
            { FragAtomId[AtomIndex] = 37; return; }

            if ((nH==0) && (dX==1) && (sAl==2))
            { FragAtomId[AtomIndex] = 38; return; }

            if ((nH==0) && (dX==1) && (sAr>=1) && (sR>=1) &&(sX==0))
            { FragAtomId[AtomIndex] = 39; return; }

            // carboxyl-like fragment R-C(=X)-X is matched also if R = H
            if (((nH==0) && (dX==1) && (sR==1) && (sX==1)) ||
                    ((nH==1) && (dX==1) && (sR==0) && (sX==1)) ||
                    ((nH==0) && (tX==1) && (sR==1)) ||
                    ((nH==0) && (dX==2)) )
            { FragAtomId[AtomIndex] = 40; return; }

            if ((nH==0) && (dX==1) && (sX==2))
            { FragAtomId[AtomIndex] = 41; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aX==1) && (asX==1))
            { FragAtomId[AtomIndex] = 42; return; }

            if ((nH==0) && (VD==3) && ((sR==1)||(aR==1)) && (aX==1) && (asX==1))
            { FragAtomId[AtomIndex] = 43; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aX==1) && (asX==1))
            { FragAtomId[AtomIndex] = 44; return; }

            ///////////////////////////////////////////////////////


            return;
        }



        //// Atom: O ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 8) {

            // Checks for particular O fragments linked to C or N

            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex)
                    continue;

                // Bound to C

                if ((ConnAugMatrix[AtomIndex][j] > 0) && (ConnAugMatrix[j][j] == 6)) {
                    boolean c_Arom = AtomAromatic[j];
                    int c_VD = 0, c_dR=0, c_dO=0, c_dX=0;
                    for (int k=0; k<nSK; k++) {
                        if (k==j)
                            continue;
                        if ( (ConnAugMatrix[j][k] > 0) && (ConnAugMatrix[k][k] !=1 ))  {
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
                    { FragAtomId[AtomIndex] = 60; return; }


                    // Particular OH groups

                    if (nH == 1) {
                        // Enol
                        if ((c_VD==3) && (!c_Arom) && (c_dR==1))
                        { FragAtomId[AtomIndex] = 57; return; }

                        // Phenol
                        if (c_Arom)
                        { FragAtomId[AtomIndex] = 57; return; }

                        // Carboxyl (matched also if R = H)
                        if ((!c_Arom) && (c_VD<=3) && (c_dO==1))
                        { FragAtomId[AtomIndex] = 57; return; }

                    }

                }


                // Bound to N

                if ((ConnAugMatrix[AtomIndex][j]>0) && (ConnAugMatrix[j][j]==7)) {
                    int n_VD = 0, n_dO=0, n_sOminus=0;
                    for (int k=0; k<nSK; k++) {
                        if (k==j)
                            continue;
                        if ( (ConnAugMatrix[j][k]>0) && (ConnAugMatrix[k][k] != 1))  {
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
                                        if ( (ConnAugMatrix[k][z] > 0) && (ConnAugMatrix[z][z] != 1)) o_VD++;
                                    }
                                    if (o_VD == 1) {
                                        int nH_O = 0;
                                        if (ExplicitHydrogen) {
                                            for (IAtom connAt : CurMol.getConnectedAtomsList(CurMol.getAtom(k)))
                                                if (connAt.getAtomicNumber() == 1)
                                                    nH_O++;
                                        } else {
                                            nH_O = CurMol.getAtom(k).getImplicitHydrogenCount();
                                        }
                                        if (nH_O == 0)
                                            n_sOminus++;
                                    }
                                }
                            }
                        }
                    }


                    // Nitro group
                    // in NO2 both oxygens are seen as O--

                    if ((VD==1) && (nH==0) && (nDouble==1) && (n_VD==3) && (n_sOminus==1))
                    { FragAtomId[AtomIndex] = 61; return; }

                    if ((VD==1) && (nH==0) && (n_VD==3) && (n_dO==1))
                    { FragAtomId[AtomIndex] = 61; return; }

                }

            }

            // correction for O-60 as in Dragon (consider also direct link to -X- in R-O-C=X)
            int nXVD2=0;
            for (int j=0; j<nSK; j++) {
                if (j == AtomIndex)
                    continue;
                if ((ConnAugMatrix[AtomIndex][j] > 0) && (IsAtomElectronegative((int)ConnAugMatrix[j][j])) ) {
                    int xVD = 0;
                    for (int k=0; k<nSK; k++) {
                        if (j == k) continue;
                        if ((ConnAugMatrix[k][j] > 0))
                            xVD++;
                    }
                    if (xVD == 2)
                        nXVD2++;
                }
            }
            if (((VD+nH)==2) && (nSingle==2) && (nXVD2==1))
            { FragAtomId[AtomIndex] = 60; return; }


            // generic alcohol (even if not bound to C)
            if ((nSingle==1) && (nH==1))
            { FragAtomId[AtomIndex] = 56; return; }


            if ((VD==1) && (nDouble==1))
            { FragAtomId[AtomIndex] = 58; return; }

            if ((VD==2) && (nSingle==2)) {
                for (int j=0; j<nSK; j++) {
                    if (j==AtomIndex)
                        continue;
                    if ((ConnAugMatrix[AtomIndex][j]==1) && (ConnAugMatrix[j][j]==8)) {
                        // Found an oxygen, R-O-O
                        FragAtomId[AtomIndex] = 63;
                        return;
                    }
                }
            }

            if (((VD==2) && (sAl==1) && (sAr==1)) ||
                    ((VD==2) && (sAr==2)) ||
                    ((VD==2) && (aR==2)))
            { FragAtomId[AtomIndex] = 60; return; }

            if ((VD==2) && (nSingle==2) && (sAr==0))
            { FragAtomId[AtomIndex] = 59; return; }

            if (Charge == -2)
            { FragAtomId[AtomIndex] = 61; return; }

            if (Charge == -1)
            { FragAtomId[AtomIndex] = 62; return; }

            return;
        }



        //// Atom: Se ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 34) {

            if ((VD==2) && (nSingle==2))
            { FragAtomId[AtomIndex] = 64; return; }

            if ((VD==1) && (nDouble==1))
            { FragAtomId[AtomIndex] = 65; return; }

            return;
        }



        //// Atom: N ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 7) {

            // Checks for particular aromatic form
            boolean AromPyridineLike=false, AromPyrroleLike=false;
            if (nArom >= 2) {
                if ((VD + nH - Charge) == 2 )
                    AromPyridineLike = true;
                if ((VD + nH - Charge) == 3 )
                    AromPyrroleLike = true;
            }

            // Checks for particular N fragments linked to C, O or X

            int nO=0, sRCO=0, sXdX=0, sXvd2dX=0;

            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex)
                    continue;
                if (ConnAugMatrix[AtomIndex][j]>0) {
                    int Z = (int)ConnAugMatrix[j][j];

                    // O
                    if (Z==8)
                        nO++;

                    // C for -RCO
                    if ( (Z==6) && (ConnAugMatrix[AtomIndex][j] == 1) )  {
                        for (int k=0; k<nSK; k++) {
                            if ((k==j) || (k==AtomIndex)) continue;
                            if ((ConnAugMatrix[k][j] == 2) && ( (ConnAugMatrix[k][k]==8) || (ConnAugMatrix[k][k]==16) ) ) {
                                sRCO++;
                                break;
                            }
                        }
                    }

                    // Electronegative
                    // NOTE: to be consistent with Dragon, also R is used
                    // instead of only electronegative atoms
                    //
                    // NOT ANYMORE (delete comment)
                    if ((IsAtomElectronegative(Z))) {
                        int x_VD=0, x_dX=0;
                        for (int k=0; k<nSK; k++) {
                            if ((k==j) || (k==AtomIndex)) continue;
                            if ((ConnAugMatrix[k][j]>0) && (ConnAugMatrix[k][k] != 1)) {
                                x_VD++;
                                if (IsAtomElectronegative((int)ConnAugMatrix[k][k]))
                                    if (ConnAugMatrix[k][j] == 2)
                                        x_dX++;
                            }
                        }
                        if ( (x_dX==1))
                            sXdX++;
                        if ( (x_dX==1) && (x_VD==2))
                            sXvd2dX++;
                    }

                }
            }

            // fragments with possible [N+](=O)[O-] or R-O-N=O
            if (((VD==3) && (nH==0) && (sAr==1) && (nO==2)) ||
//                    ((VD==3) && (nH==0) && (aR==2) && (aX==1)) ||
                    ((VD==3) && (nH==0) && (aR==2) && (nO==1)) ||
                    ((VD==2) && (nH==0) && (nO==2)))
//                || ((VD==2) && (nH==0) && (nO==2)))
            { FragAtomId[AtomIndex] = 76; return; }

            if ((VD==3) && (nH==0) && (sAr==0) && (nO>=2))
            { FragAtomId[AtomIndex] = 77; return; }

            // N Charged +1
            if ((!AromPyrroleLike) && (!AromPyridineLike))
                if (Charge==1)
                { FragAtomId[AtomIndex] = 79; return; }

            // fragment with particular groups
            if (!AromPyrroleLike)
                if ((((VD+nH)==3) && (sRCO>0)) ||   // RCO-N<
                        (((VD+nH)==3) && (sXvd2dX>0)))     // >N-X=X
                { FragAtomId[AtomIndex] = 72; return; }


            if ((VD==1) && (nH==2) && (sAr==0) && (sX==0) && (nSingle==1))
            { FragAtomId[AtomIndex] = 66; return; }

            if ((VD==2) && (nH==1) && (sAr==0) && (nSingle==2))
            { FragAtomId[AtomIndex] = 67; return; }

            if ((VD==3) && (nH==0) && (sAr==0) && (aAr==0) && (nSingle==3))
            { FragAtomId[AtomIndex] = 68; return; }

            if ((VD==1) && (nH==2) && ( (sAr==1) || (sX==1) ))
            { FragAtomId[AtomIndex] = 69; return; }

            if ((VD==2) && (nH==1) && (sAl==1) && (sAr==1))
            { FragAtomId[AtomIndex] = 70; return; }

            if ((VD==3) && (nH==0) && (sAl==2) && (sAr==1))
            { FragAtomId[AtomIndex] = 71; return; }

            if (((VD==2) && (nH==1) && (sAr==2)) ||
                    ((VD==3) && (nH==0) && (sAr==3)) ||
                    ((VD==3) && (nH==0) && (sAr==2) && (sAl==1)) ||
                    ((AromPyrroleLike)))
            { FragAtomId[AtomIndex] = 73; return; }

            // R=N- is matched also if it is R=NH
            if (((VD==1) && (nH==0) && (tR==1)) ||
                    ((VD==2) && (nH==0) && (dR==1) && (nSingle==1)) ||
                    ((VD==1) && (nH==1) && (dR==1)))
            { FragAtomId[AtomIndex] = 74; return; }

            if (((AromPyridineLike)) ) // ||
//                ((VD==2) && (aR==1) && (aX==1)))
            { FragAtomId[AtomIndex] = 75; return; }

            if (((VD==2) && (nH==0) && (sAr==1) && (dX==1)) ||
                    ((VD==2) && (nH==0) && (sX==1) && (dX==1)))
            { FragAtomId[AtomIndex] = 78; return; }

            return;
        }



        //// Halogen ions //////////////////////////////////////////////////////

        // F
        if ((ConnAugMatrix[AtomIndex][AtomIndex] == 9) && (VD==0))
        { FragAtomId[AtomIndex] = 101; return; }

        // Cl
        if ((ConnAugMatrix[AtomIndex][AtomIndex] == 17) && (VD==0))
        { FragAtomId[AtomIndex] = 102; return; }

        // Br
        if ((ConnAugMatrix[AtomIndex][AtomIndex] == 35) && (VD==0))
        { FragAtomId[AtomIndex] = 103; return; }

        // I
        if ((ConnAugMatrix[AtomIndex][AtomIndex] == 53) && (VD==0))
        { FragAtomId[AtomIndex] = 104; return; }



        //// Atom: S ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 16) {

            int nO=0;

            for (int j=0; j<nSK; j++) {
                if (j==AtomIndex) continue;
                if (ConnAugMatrix[AtomIndex][j]>0) {
                    if (ConnAugMatrix[j][j]==8)
                        nO++;
                }
            }

            if ((nH==1) && (sR==1))
            { FragAtomId[AtomIndex] = 106; return; }

            if ((VD==2) && (nH==0))    // RSR is valid both for single and aromatic bonds (?)
            { FragAtomId[AtomIndex] = 107; return; }

            if ((VD==1) && (nDouble==1))
            { FragAtomId[AtomIndex] = 108; return; }

            if ((VD==3) && (nSingle==2) && (nO==1))
            { FragAtomId[AtomIndex] = 109; return; }

            if ((VD==4) && (nSingle==2) && (nO>=2))
            { FragAtomId[AtomIndex] = 110; return; }

            return;
        }



        //// Atom: Si //////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 14) {

            if ((VD==4) && (nSingle==4))
            { FragAtomId[AtomIndex] = 111; return; }

            return;

        }



        //// Atom: B //////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 5) {

            if ((VD==3) && (nSingle==3))
            { FragAtomId[AtomIndex] = 112; return; }

            return;

        }



        //// Atom: P //////////////////////////////////////////////////////////

        if (ConnAugMatrix[AtomIndex][AtomIndex] == 15) {

            if ((VD==4) && (Charge==1))
            { FragAtomId[AtomIndex] = 115; return; }

            if ((VD==4) && (dX==1) && (sR==3))
            { FragAtomId[AtomIndex] = 116; return; }

            if ((VD==4) && (dX==1) && (sX==3))
            { FragAtomId[AtomIndex] = 117; return; }


            if ((VD==3) && (sX==3))
            { FragAtomId[AtomIndex] = 118; return; }

            if ((VD==3) && (sR==3))
            { FragAtomId[AtomIndex] = 119; return; }

            // consider as phosponate also PH(=X)(X)(X)
            if ( ((VD==4) && (dX==1) && (sR==1) && (sX==2)) ||
                    ((VD==3) && (dX==1) && (sX==2)))
            { FragAtomId[AtomIndex] = 120; }

        }
    }


    /**
     * Return true if an atom is electronegative, i.e. one of the following:
     * O, N, S, P, B, Si, Se or Halogens (F, Cl, Br, I)
     *
     * @param AtomicNumber
     * @return
     */
    private static boolean IsAtomElectronegative(int AtomicNumber) {

        // O, N, S, P, B, Si, Se, halogen

        if ((AtomicNumber==7)||(AtomicNumber==8)||(AtomicNumber==15)||
                (AtomicNumber==16)||(AtomicNumber==34)||(AtomicNumber==9)||
                (AtomicNumber==5)||(AtomicNumber==14)||
                (AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53)) {
            return true;
        }
        return false;
    }

    /**
     * Return true if an atom is Halogen (F, Cl, Br, I)
     *
     * @param AtomicNumber
     * @return
     */
    private static boolean IsAtomHalogen(int AtomicNumber) {
        if ((AtomicNumber==9)||(AtomicNumber==17)||(AtomicNumber==35)||(AtomicNumber==53))
            return true;
        return false;
    }

}
