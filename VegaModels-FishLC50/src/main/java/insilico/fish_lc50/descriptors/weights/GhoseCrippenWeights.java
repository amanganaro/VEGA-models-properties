package insilico.fish_lc50.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

public class GhoseCrippenWeights {

    final static double[] Hydrophobicity = {
            0,  // foo for index 0
            -1.5603,
            -1.012,
            -0.6681,
            -0.3698,
            -1.788,
            -1.2486,
            -1.0305,
            -0.6805,
            -0.3858,
            0.7555,
            -0.2849,
            0.02,
            0.7894,
            1.6422,
            -0.7866,
            -0.3962,
            0.0383,
            -0.8051,
            -0.2129,
            0.2432,
            0.4697,
            0.2952,
            0,
            -0.3251,
            0.1492,
            0.1539,
            0.0005,
            0.2361,
            0.3514,
            0.1814,
            0.0901,
            0.5142,
            -0.3723,
            0.2813,
            0.1191,
            -0.132,
            -0.0244,
            -0.2405,
            -0.0909,
            -0.1002,
            0.4182,
            -0.2147,
            -0.0009,
            0.1388,
            0,
            0.7341,
            0.6301,
            0.518,
            -0.0371,
            -0.1036,
            0.5234,
            0.6666,
            0.5372,
            0.6338,
            0.362,
            -0.3567,
            -0.0127,
            -0.0233,
            -0.1541,
            0.0324,
            1.052,
            -0.7941,
            0.4165,
            0.6601,
            0,
            -0.5427,
            -0.3168,
            0.0132,
            -0.3883,
            -0.0389,
            0.1087,
            -0.5113,
            0.1259,
            0.1349,
            -0.1624,
            -2.0585,
            -1.915,
            0.4208,
            -1.4439,
            0,
            0.4797,
            0.2358,
            0.1029,
            0.3566,
            0.1988,
            0.7443,
            0.5337,
            0.2996,
            0.8155,
            0.4856,
            0.8888,
            0.7452,
            0.5034,
            0.8995,
            0.5946,
            1.4201,
            1.1472,
            0,
            0.7293,
            0.7173,
            0,
            -2.6737,
            -2.4178,
            -3.1121,
            0,
            0.6146,
            0.5906,
            0.8758,
            -0.4979,
            -0.3786,
            1.5188,
            1.0255,
            0,
            0,
            0,
            -0.9359,
            -0.1726,
            -0.7966,
            0.6705,
            -0.4801 };

    final static double[] MolarRefractivity = {
            0,  // foo for index 0
            2.968,
            2.9116,
            2.8028,
            2.6205,
            3.015,
            2.9244,
            2.6329,
            2.504,
            2.377,
            2.5559,
            2.303,
            2.3006,
            2.9627,
            2.3038,
            3.2001,
            4.2654,
            3.9392,
            3.6005,
            4.487,
            3.2001,
            3.4825,
            4.2817,
            3.9556,
            3.4491,
            3.8821,
            3.7593,
            2.5009,
            2.5,
            3.0627,
            2.5009,
            0,
            2.6632,
            3.4671,
            3.6842,
            2.9372,
            4.019,
            4.777,
            3.9031,
            3.9964,
            3.4986,
            3.4997,
            2.7784,
            2.6267,
            2.5,
            0,
            0.8447,
            0.8939,
            0.8005,
            0.832,
            0.8,
            0.8188,
            0.9215,
            0.9769,
            0.7701,
            0,
            1.7646,
            1.4778,
            1.4429,
            1.6191,
            1.3502,
            1.945,
            0,
            0,
            11.1366,
            13.1149,
            2.6221,
            2.5,
            2.898,
            3.6841,
            4.2808,
            3.6189,
            2.5,
            2.7956,
            2.7,
            4.2063,
            4.0184,
            3.0009,
            4.7142,
            0,
            0,
            0.8725,
            1.1837,
            1.1573,
            0.8001,
            1.5013,
            5.6156,
            6.1022,
            5.9921,
            5.3885,
            6.1363,
            8.5991,
            8.9188,
            8.8006,
            8.2065,
            8.7352,
            13.9462,
            14.0792,
            14.073,
            12.9918,
            13.3408,
            0,
            0,
            0,
            0,
            0,
            7.8916,
            7.7935,
            9.4338,
            7.7223,
            5.7558,
            0,
            0,
            0,
            0,
            0,
            5.5306,
            5.5152,
            6.836,
            10.0101,
            5.2806 };


    /**
     * Provides the hydrophobicity contribution for an atom centered fragment
     *
     * @param FragmentIndex index of the ACF
     * @return hydrophobicity contribution
     */
    public static double GetHydrophobiticty(int FragmentIndex) {
        if ((FragmentIndex>=1) && (FragmentIndex<=120))
            return Hydrophobicity[FragmentIndex];
        else
            return 0;

    }

    /**
     * Provides the molar refractivity contribution for an atom centered fragment
     *
     * @param FragmentIndex index of the ACF
     * @return molar regfractivity contribution
     */
    public static double GetMolarRefractivity(int FragmentIndex) {
        if ((FragmentIndex>=1) && (FragmentIndex<=120))
            return MolarRefractivity[FragmentIndex];
        else
            return 0;
    }

    private static boolean IsAtomElectronegative(int AtomicNumber) {
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

    public static double[] GetHydrophobiticty(IAtomContainer CurMol) {

        int[] FragAtomId = new int[CurMol.getAtomCount()];
        double[] w = new double[CurMol.getAtomCount()];
        for (int i=0; i<CurMol.getAtomCount(); i++) {
            FragAtomId[i] = Descriptor.MISSING_VALUE;
            w[i] = Descriptor.MISSING_VALUE;
        }

        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
        } catch (Exception e) {
            return w;
        }

        for (IAtom at : CurMol.atoms())
            GetHydrophobitictyPrivate(at, CurMol, FragAtomId, ConnAugMatrix);

        for (int i=0; i<CurMol.getAtomCount(); i++) {
            if (FragAtomId[i] != Descriptor.MISSING_VALUE)
                w[i] = Hydrophobicity[FragAtomId[i]];
        }

        return w;
    }

    public static double[] GetMolarRefractivity(IAtomContainer CurMol) {

        int[] FragAtomId = new int[CurMol.getAtomCount()];
        double[] w = new double[CurMol.getAtomCount()];
        for (int i=0; i<CurMol.getAtomCount(); i++) {
            FragAtomId[i] = Descriptor.MISSING_VALUE;
            w[i] = Descriptor.MISSING_VALUE;
        }

        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = ConnectionAugMatrix.getMatrix(CurMol);
        } catch (Exception e) {
            return w;
        }

        for (IAtom at : CurMol.atoms())
            GetHydrophobitictyPrivate(at, CurMol, FragAtomId, ConnAugMatrix);

        for (int i=0; i<CurMol.getAtomCount(); i++) {
            if (FragAtomId[i] != Descriptor.MISSING_VALUE)
                w[i] = MolarRefractivity[FragAtomId[i]];
        }

        return w;
    }


    private static void GetHydrophobitictyPrivate(IAtom at, IAtomContainer CurMol, int[] FragAtomId, double[][] ConnAugMatrix) {

        int At = CurMol.indexOf(at);
        int nSK =  CurMol.getAtomCount();

        boolean[] AtomAromatic = new boolean[nSK];
        for (int i=0; i<nSK; i++)
            if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC))
                AtomAromatic[i] = true;
            else
                AtomAromatic[i] = false;


        ////////////////////////////
        // taken from ACF code
        ////////////////////////////

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

            C_OxiNumber += (sX * 1) + (dX * 2) + (tX * 3);
            C_OxiNumber += (asX * 1);
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

                    H_type=nH;

                } else {

                    // C0sp3 (no X attached to next C)
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==0))
                        H_type=nH;

                    // C1sp3, C0sp2
                    if ( ((C_OxiNumber==1) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==2)) )
                        H_type=nH;

                    // C2sp3, C1sp2, C0sp
                    if ( ((C_OxiNumber==2) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==1) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==0) && (C_Hybridazion==1)) )
                        H_type=nH;

                    // C3sp3, C2sp2, C2sp2, C3sp
                    if ( ((C_OxiNumber==3) && (C_Hybridazion==3)) ||
                            ((C_OxiNumber==2) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==2)) ||
                            ((C_OxiNumber==3) && (C_Hybridazion==1)) )
                        H_type=nH;

                    // C0sp3 with 1 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==1))
                        H_type=nH;

                    // C0sp3 with 2 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==2))
                        H_type=nH;

                    // C0sp3 with 3 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==3))
                        H_type=nH;

                    // C0sp3 with 4 X atom attached to next C
                    if ((C_OxiNumber==0) && (C_Hybridazion==3) && (C_CX==4))
                        H_type=nH;

                }

            } else {

                // H to heteroatom
                H_type=nH;

            }

            for (IAtom conn : CurMol.getConnectedAtomsList(at))
                if (conn.getSymbol().equalsIgnoreCase("H"))
                    FragAtomId[CurMol.indexOf(conn)] = H_type;

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

                    if (ConnAugMatrix[At][At] == 6) {
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

                    if (ConnAugMatrix[At][At] == 6) {
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

                    if (ConnAugMatrix[At][At] == 6) {
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

        if (ConnAugMatrix[At][At] == 6) {

            if ( ((nH==3) && (sR==1)) || (nH==4) )
            { FragAtomId[At] = 1; return; }

            if ((nH==2) && (sR==2))
            { FragAtomId[At] = 2; return; }

            if ((nH==1) && (sR==3))
            { FragAtomId[At] = 3; return; }

            if ((nH==0) && (sR==4))
            { FragAtomId[At] = 4; return; }

            if ((nH==3) && (sX==1))
            { FragAtomId[At] = 5; return; }

            if ((nH==2) && (sR==1) && (sX==1))
            {FragAtomId[At] = 6; return; }

            if ((nH==2) && (sX==2))
            { FragAtomId[At] = 7; return; }

            if ((nH==1) && (sR==2) && (sX==1))
            { FragAtomId[At] = 8; return; }

            if ((nH==1) && (sR==1) && (sX==2))
            { FragAtomId[At] = 9; return; }

            if ((nH==1) && (sX==3))
            { FragAtomId[At] = 10; return; }

            if ((nH==0) && (sR==3) && (sX==1))
            { FragAtomId[At] = 11; return; }

            if ((nH==0) && (sR==2) && (sX==2))
            { FragAtomId[At] = 12; return; }

            if ((nH==0) && (sR==1) && (sX==3))
            { FragAtomId[At] = 13; return; }

            if ((nH==0) && (sX==4))
            { FragAtomId[At] = 14; return; }

            if ((nH==2) && (dR==1))
            { FragAtomId[At] = 15; return; }

            if ((nH==1) && (dR==1) && (sR==1))
            { FragAtomId[At] = 16; return; }

            if ((nH==0) && (dR==1) && (sR==2))
            { FragAtomId[At] = 17; return; }

            if ((nH==1) && (dR==1) && (sX==1))
            { FragAtomId[At] = 18; return; }

            if ((nH==0) && (dR==1) && (sX==1) && (sR==1))
            { FragAtomId[At] = 19; return; }

            if ((nH==0) && (dR==1) && (sX==2))
            { FragAtomId[At] = 20; return; }

            if ((nH==1) && (tR==1))
            { FragAtomId[At] = 21; return; }

            if ( ((nH==0) && (tR==1) && (sR==1)) ||
                    ((nH==0) && (dR==2)) )
            { FragAtomId[At] = 22; return; }

            if ((nH==0) && (tR==1) && (sX==1))
            { FragAtomId[At] = 23; return; }

            if ((nH==1) && (VD==2) && (aR==2))
            { FragAtomId[At] = 24; return; }

            if ((nH==0) && (VD==3) && (aR>=2) && ((sR==1)||(aR==3)))
            { FragAtomId[At] = 25; return; }

            if ((nH==0) && (VD==3) && (aR==2) && (sX==1))
            { FragAtomId[At] = 26; return; }

            if ((nH==1) && (VD==2) && (aR==1) && (aX==1))
            { FragAtomId[At] = 27; return; }

            if ((nH==0) && (VD==3) && (aR>=1) && (aX==1) && ((sR==1)||(aR==2)))
            { FragAtomId[At] = 28; return; }

            if ((nH==0) && (VD==3) && (aR==1) && (aX==1) && (sX==1))
            { FragAtomId[At] = 29; return; }

            if ((nH==1) && (VD==2) && (aX==2))
            { FragAtomId[At] = 30; return; }

            if ((nH==0) && (VD==3) && (aX==2) && ((sR==1)||(aR==1)))
            { FragAtomId[At] = 31; return; }

            if ((nH==0) && (VD==3) && (aX==2) && (sX==1))
            { FragAtomId[At] = 32; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aR==1) && (asX==1))
            { FragAtomId[At] = 33; return; }

            if ((nH==0) && (VD==3) && (((sR==1) && (aR==1))||(aR==2)) && (asX==1))
            { FragAtomId[At] = 34; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aR==1) && (asX==1))
            { FragAtomId[At] = 35; return; }

            ///////////////////////////////////////////////////////


            if ((nH==1) && (dX==1) && (sAl==1))
            { FragAtomId[At] = 36; return; }

            if ((nH==1) && (dX==1) && (sAr==1))
            { FragAtomId[At] = 37; return; }

            if ((nH==0) && (dX==1) && (sAl==2))
            { FragAtomId[At] = 38; return; }

            if ((nH==0) && (dX==1) && (sAr==1) && (sR>=1) &&(sX==0))
            { FragAtomId[At] = 39; return; }

            if (((nH==0) && (dX==1) && (sR==1) && (sX==1)) ||
                    ((nH==0) && (tX==1) && (sR==1)) ||
                    ((nH==0) && (dX==2)) )
            { FragAtomId[At] = 40; return; }

            if ((nH==0) && (dX==1) && (sX==2))
            { FragAtomId[At] = 41; return; }


            ///// da controllare (differenza fra -- e ..) ////////

            if ((nH==1) && (VD==2) && (aX==1) && (asX==1))
            { FragAtomId[At] = 42; return; }

            if ((nH==0) && (VD==3) && ((sR==1)||(aR==1)) && (aX==1) && (asX==1))
            { FragAtomId[At] = 43; return; }

            if ((nH==0) && (VD==3) && (sX==1) && (aX==1) && (asX==1))
            { FragAtomId[At] = 44; return; }

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
                    { FragAtomId[At] = 60; return; }


                    // Particular OH groups

                    if (nH == 1) {
                        // Enol
                        if ((c_VD==3) && (!c_Arom) && (c_dR==1))
                        { FragAtomId[At] = 57; return; }

                        // Phenol
                        if (c_Arom)
                        { FragAtomId[At] = 57; return; }

                        // Carboxyl
                        if ((!c_Arom) && (c_VD==3) && (c_dO==1))
                        { FragAtomId[At] = 57; return; }

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
                    { FragAtomId[At] = 61; return; }

                    if ((VD==1) && (nH==0) && (n_VD==3) && (n_dO==1))
                    { FragAtomId[At] = 61; return; }

                }

            }

            // generic alcohol (even if not bound to C)
            if ((nSingle==1) && (nH==1))
            { FragAtomId[At] = 56; return; }


            if ((VD==1) && (nDouble==1))
            { FragAtomId[At] = 58; return; }

            if ((VD==2) && (nSingle==2)) {
                for (int j=0; j<nSK; j++) {
                    if (j==At)
                        continue;
                    if ((ConnAugMatrix[At][j]==1) && (ConnAugMatrix[j][j]==8)) {
                        // Found an oxygen, R-O-O
                        FragAtomId[At] = 63;
                        return;
                    }
                }
            }

            if (((VD==2) && (sAl==1) && (sAr==1)) ||
                    ((VD==2) && (sAr==2)) ||
                    ((VD==2) && (aR==2)))
            { FragAtomId[At] = 60; return; }

            if ((VD==2) && (nSingle==2) && (sAr==0))
            { FragAtomId[At] = 59; return; }

            if (Charge == -2)
            { FragAtomId[At] = 61; return; }

            if (Charge == -1)
            { FragAtomId[At] = 62; return; }

            return;
        }



        //// Atom: Se ///////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 34) {

            if ((VD==2) && (nSingle==2))
            { FragAtomId[At] = 64; return; }

            if ((VD==1) && (nDouble==1))
            { FragAtomId[At] = 65; return; }

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
            { FragAtomId[At] = 76; return; }

            if ((VD==3) && (nH==0) && (sAl==1) && (nO==2))
            { FragAtomId[At] = 77; return; }

            // N Charged +1
            if ((!AromPyrroleLike) && (!AromPyridineLike))
                if (Charge==1)
                { FragAtomId[At] = 79; return; }

            // fragment with particular groups
            if (!AromPyrroleLike)
                if (((VD==3) && (sOCR==1)) ||   // >N-OCR
                        (((VD+nH)==3) && (sXdX>0)))     // >N-X=X  NOTE: also >N-R=X for Dragon compatibility
                { FragAtomId[At] = 72; return; }


            if ((VD==1) && (nH==2) && (sAl==1))
            { FragAtomId[At] = 66; return; }

            if ((VD==2) && (nH==1) && (sAl==2))
            { FragAtomId[At] = 67; return; }

            if ((VD==3) && (nH==0) && (sAl==3))
            { FragAtomId[At] = 68; return; }

            if ((VD==1) && (nH==2) && ( (sAr==1) || (sX==1) ))
            { FragAtomId[At] = 69; return; }

            if ((VD==2) && (nH==1) && (sAl==2))
            { FragAtomId[At] = 70; return; }

            if ((VD==3) && (nH==0) && (sAl==2) && (sAr==1))
            { FragAtomId[At] = 71; return; }

            if (((VD==2) && (nH==1) && (sAr==2)) ||
                    ((VD==3) && (nH==0) && (sAr==3)) ||
                    ((VD==3) && (nH==0) && (sAr==2) && (sAl==1)) ||
                    ((AromPyrroleLike)))
            { FragAtomId[At] = 73; return; }

            if (((VD==1) && (nH==0) && (tR==1)) ||
                    ((VD==2) && (nH==0) && (dR==1) && (sR==1)))
            { FragAtomId[At] = 74; return; }

            if (((AromPyridineLike)) ) // ||
//                ((VD==2) && (aR==1) && (aX==1)))
            { FragAtomId[At] = 75; return; }

            if (((VD==2) && (nH==0) && (sAr==1) && (dX==1)) ||
                    ((VD==2) && (nH==0) && (sX==1) && (dX==1)))
            { FragAtomId[At] = 78; return; }

            return;
        }



        //// Halogen ions //////////////////////////////////////////////////////

        // F
        if ((ConnAugMatrix[At][At] == 9) && (VD==0))
        { FragAtomId[At] = 101; return; }

        // Cl
        if ((ConnAugMatrix[At][At] == 17) && (VD==0))
        { FragAtomId[At] = 102; return; }

        // Br
        if ((ConnAugMatrix[At][At] == 35) && (VD==0))
        { FragAtomId[At] = 103; return; }

        // I
        if ((ConnAugMatrix[At][At] == 53) && (VD==0))
        { FragAtomId[At] = 104; return; }



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
            { FragAtomId[At] = 106; return; }

            if ((VD==2) && (nH==0))    // RSR is valid both for single and aromatic bonds (?)
            { FragAtomId[At] = 107; return; }

            if ((VD==1) && (nDouble==1))
            { FragAtomId[At] = 108; return; }

            if ((VD==3) && (nSingle==2) && (nO==1))
            { FragAtomId[At] = 109; return; }

            if ((VD==4) && (nSingle==2) && (nO>=2))
            { FragAtomId[At] = 110; return; }

            return;
        }



        //// Atom: Si //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 14) {

            if ((VD==4) && (nSingle==4))
            { FragAtomId[At] = 111; return; }

            return;

        }



        //// Atom: B //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 5) {

            if ((VD==3) && (nSingle==3))
            { FragAtomId[At] = 112; return; }

            return;

        }



        //// Atom: P //////////////////////////////////////////////////////////

        if (ConnAugMatrix[At][At] == 15) {

            if ((VD==4) && (Charge==1))
            { FragAtomId[At] = 115; return; }

            if ((VD==4) && (dX==1) && (sR==3))
            { FragAtomId[At] = 116; return; }

            if ((VD==4) && (dX==1) && (sX==3))
            { FragAtomId[At] = 117; return; }


            if ((VD==3) && (sX==3))
            { FragAtomId[At] = 118; return; }

            if ((VD==3) && (sR==3))
            { FragAtomId[At] = 119; return; }

            if ((VD==4) && (dX==1) && (sR==1) && (sX==2))
            { FragAtomId[At] = 120; return; }

        }


    }
}
