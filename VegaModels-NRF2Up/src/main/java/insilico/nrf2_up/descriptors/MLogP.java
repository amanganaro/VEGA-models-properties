package insilico.nrf2_up.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.InsilicoMoleculeNormalization;
import insilico.core.tools.utils.MoleculeUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;

import java.util.Arrays;
import java.util.List;

/**
 * MLogP (Moriguchi LogP) descriptor.
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class MLogP extends DescriptorBlock {
    private static final Logger log = LogManager.getLogger(MLogP.class);

    private final static long serialVersionUID = 1L;

    private final static String BlockName = "MLOGP";

    private IAtomContainer CurMol;
    private int nAtoms;
    private double[][] ConnMatrix;
    private int[][] TopoMatrix;
    private RingSet MolRings;


    /**
     * Constructor.
     */
    public MLogP() {
        super();
        this.Name = MLogP.BlockName;
    }


    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        this.Add("MLogP","");
        SetAllValues(Descriptor.MISSING_VALUE);
    }


    /**
     * Calculate descriptors for the given molecule.
     *
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {

        GenerateDescriptors();

        try {
            CurMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        // Gets matrices
        try {
            ConnMatrix = mol.GetMatrixConnectionAugmented();
            TopoMatrix = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        // !!! in origine usava topo matrix di CDK

        double CX = 0;
        int NO = 0;
        int PRX = 0;
        int UB = 0;
        int HB = 0;
        int POL = 0;
        double AMP = 0;
        int ALK = 0;
        int RNG = 0;
        double QN = 0;
        int NO2 = 0;
        double NCS = 0;
        int BLM = 0;

        double LogP = 0;

        nAtoms = CurMol.getAtomCount();
        try {
            MolRings = mol.GetSSSR();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }


        // Atom Counters

        int nC=0, nN=0, nO=0, nF=0, nP=0, nS=0, nCl=0, nBr=0, nI=0;
        int nNO2=0, nNamm=0, nNoxide=0;
        for (int i=0; i<nAtoms; i++) {
            int atom = (int)ConnMatrix[i][i];
            switch (atom) {
                case 6: nC++; break;
                case 7: nN++; break;
                case 8: nO++; break;
                case 9: nF++; break;
                case 15: nP++; break;
                case 16: nS++; break;
                case 17: nCl++; break;
                case 35: nBr++; break;
                case 53: nI++; break;
            }

            // N centered groups
            if (atom == 7) {

                int nOdbl = 0;
                int nOminus = 0;
                int nVD = 0;
                int nCoxide = 0;
                double totBondOrder = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (j == i) continue;
                    if (ConnMatrix[i][j] > 0) {
                        nVD++;
                        totBondOrder += ConnMatrix[i][j];
                        if  (ConnMatrix[j][j] == 6)
                            nCoxide++;
                        if ((ConnMatrix[i][j] == 2) && (ConnMatrix[j][j] == 8))
                            nOdbl++;
                        if ((ConnMatrix[i][j] == 1) && (ConnMatrix[j][j] == 8)) {
                            int OVD = 0;
                            for (int k = 0; k < nAtoms; k++) {
                                if (k == j) continue;
                                if (ConnMatrix[j][k] > 0) OVD++;
                            }
                            if (OVD == 1)
                                nOminus++;
                        }
                    }
                }
                int nH;
                try {
                    nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }

                // count quaternary N (or N oxide with VD = 4)
                if ((nVD+nH) == 4) {
                    if ((nOdbl == 1) && (nCoxide + nH == 3))
                        nNoxide++;
                    else
                        nNamm++;
                }

                // count N oxide with VD = 2
                if ( (nVD == 2) && (totBondOrder == 4) )
                    if (nOdbl > 0)
                        nNoxide++;

                // N oxide in pyridine type structure
                if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC))
                    if ( (nVD == 3) && (nOdbl == 1) && (nCoxide == 2) )
                        nNoxide++;

                // in Dragon seems to ignore other possible N oxides ? possible error

                // count NO2 groups
                if (InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION) {
                    // due to preprocessing, nitro groups are  set in O=N=O form
                    if ((nVD == 3) && (nOdbl == 2) )
                        nNO2++;
                } else {
                    // due to preprocessing, nitro groups are  set in O=[N+][O-] form
                    if ((nVD == 3) && (nOdbl == 1) && (nOminus == 1))
                        nNO2++;
                }
            }
        }


        //// 1 - CX

        CX = nC + 0.5 * nF + nCl + 1.5 * nBr + 2 * nI;


        //// 2 - NO

        NO = nN + nO;


        //// 3 - PRX

        PRX = Calculate_PRX();


        //// 4 - UB

        for (IBond bnd : CurMol.bonds()) {
            // new CDK correctly parses aromatic bonds to kekule structure
            // indeed their order is set correctly to single or double
            if ( (bnd.getOrder().numeric() == 2) || (bnd.getOrder().numeric() == 3) )
                UB++;
        }

        // don't count double bonds in NO2 (consider one dbl for each NO2 group with standard normalization)
        int DblBondsInNO2 = 1;
        if (InsilicoMoleculeNormalization.DRAGON7_COMPLIANT_NORMALIZATION)
            DblBondsInNO2 = 2;
        UB -= nNO2 * DblBondsInNO2;

        // Correction for aromatic N like in pyridine oxide
        for (int i=0; i<nAtoms; i++) {
            if ( (ConnMatrix[i][i] == 7) && (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) ) {
                // TODO: ma cosa fa in dragon??
            }
        }


        //// 5. HB
        HB = Check_HB() ? 1 : 0;


        //// 6. POL
        POL = Calculate_POL();
        if (POL > 4)
            POL = 4;


        //// 7. AMP
        AMP = Calculate_AMP(CurMol);


        //// 8. ALK
        ALK = Calculate_ALK(CurMol);


        //// 9. RNG
        RNG = Calculate_RNG(CurMol);


        //// 10. QN
        QN = nNamm + 0.5 * nNoxide;


        //// 11. NO2
        NO2 = nNO2;


        //// 12. NCS
        NCS = Check_NCS();


        //// 13. BLM

        for (int i=0; i<nAtoms; i++) {

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("N")) {

                if (BLM == 0)
                    BLM = Check_BLM(CurMol, i);

            }

        }



        // Calculation of LogP
        LogP = -1.041 + 1.244 * Math.pow(CX,0.6) - 1.017 * Math.pow(NO, 0.9) +
                0.406 * PRX + -0.145 * Math.pow(UB, 0.8) + 0.511 * HB +
                0.268 * POL - 2.215 * AMP + 0.912 * ALK -0.392 * RNG -
                3.684 * QN + 0.474 * NO2 + 1.582 * NCS + 0.773 * BLM;

        // Sets final descriptor value
        this.SetByName("MLogP", LogP);

    }


    private int Calculate_PRX() {

        int XAYbond = 0;

        for (int i=0; i<nAtoms; i++) {

            int nNbond = 0, nObond = 0, nOdbl = 0, nN_VD =0;

            String sym = CurMol.getAtom(i).getSymbol();

            if ( (sym.equalsIgnoreCase("C")) || (sym.equalsIgnoreCase("S")) ||
                    (sym.equalsIgnoreCase("P")) ) {

                for (int j=0; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j] != 0) {

                        if (ConnMatrix[j][j] == 7) {
                            // N
                            nNbond++;
                            for (int k=0; k<nAtoms; k++) {
                                if (k==j) continue;
                                if (ConnMatrix[j][k] != 0) nN_VD++;
                            }

                            // counts also H in N vd?
                            int nH;
                            try {
                                nH = CurMol.getAtom(j).getImplicitHydrogenCount();
                            } catch (Exception e) {
                                nH = 0;
                            }
                            nN_VD += nH;

                        }

                        if (ConnMatrix[j][j] == 8) {
                            // O
                            nObond++;
                            if (ConnMatrix[i][j] == 2) nOdbl++;
                        }
                    }
                }

                int nNObond = nObond + nNbond;

                // X-A-Y
                if ( nNObond == 2) {
                    XAYbond += 2;

                    // Corrections

                    // O-C-O
                    if ( (nObond == 2) && (nOdbl == 0) ) XAYbond--;

                    // >N-C-N<
                    if ( (nNbond == 2) && (nN_VD == 6) ) XAYbond--;

                    // >N-C-O or >N-C=O
                    if ( (nNbond == 1) && (nN_VD == 3) && (nObond == 1) ) XAYbond--;

                }

                if ( nNObond == 3) {
                    XAYbond += 3;
                    // Correction for Sulfonamide S(=O)(=O)N and Carboxamide
                    if ( (nOdbl > 0) && (((nNbond == 1) && (nN_VD == 3)) || ((nNbond == 2) && (nN_VD > 4))) ) XAYbond--;
                }

                if (nNObond == 4) {
                    XAYbond += 6;
                    if ( (sym.equalsIgnoreCase("S")) && (nOdbl > 0) &&
                            (((nNbond == 1) && (nN_VD == 3)) || ((nNbond == 2) && (nN_VD > 4)) || ((nNbond == 3)
                            && (nN_VD > 6))) )
                        XAYbond--;
                }
            }
        }

        // NO
        int nNObond = 0;
        for (int i=0; i<nAtoms-1; i++) {
            if ( (ConnMatrix[i][i] == 7) || (ConnMatrix[i][i] == 8) )
                for (int j = i+1; j < nAtoms; j++)
                    if ( (ConnMatrix[j][j] == 7) || (ConnMatrix[j][j] == 8) ) {
                        if (ConnMatrix[i][j] != 0)
                            nNObond++;
                    }
        }

        return nNObond * 2 + XAYbond;
    }


    private boolean Check_HB() {

        for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
            IRing curRing = (IRing) MolRings.getAtomContainer(r);

            // consider aromatic PSA rings as in Dragon
            // Dragon seems NOT to consider if the ring structure is aromatic
            boolean AromaticPSA = true;
            for (int i = 0; i < curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));
                int VD = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (atIdx == j) continue;
                    if (ConnMatrix[atIdx][j] != 0)
                        VD++;
                }
                int nH;
                try {
                    nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }

                VD += nH;
                if (VD > 3) {
                    AromaticPSA = false;
                    break;
                }
            }

            if (!AromaticPSA)
                continue;

            // check all ring atoms for substituent
            // to find NH2 - OH in orto position
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));

                int VD = 0;
                for (int j = 0; j < nAtoms; j++) {
                    if (atIdx == j) continue;
                    if (ConnMatrix[atIdx][j] != 0)
                        VD++;
                }
                int nH;
                try {
                    nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
                } catch (Exception e) {
                    nH = 0;
                }
                VD += nH;

                if (VD == 3) {

                    int FoundNH2 = Find_Attached_NH2(atIdx);
                    int FoundOH = Find_Attached_OH(atIdx);

                    // found a NH2 or OH group on atIdx
                    if ( (FoundOH > -1) || (FoundNH2 > -1) ) {

                        for (int atIdxB = 0; atIdxB < nAtoms; atIdxB++) {
                            if (atIdx == atIdxB) continue;

                            // check connected atoms in the same ring
                            if ( (ConnMatrix[atIdx][atIdxB] != 0) && (curRing.contains(CurMol.getAtom(atIdxB))) )  {

                                int FoundNH2_b = Find_Attached_NH2(atIdxB);
                                int FoundOH_b = Find_Attached_OH(atIdxB);

                                // check C=O not belonging to the same cycle
                                int FoundCO = -1;
                                for (int atIdxC = 0; atIdxC < nAtoms; atIdxC++) {
                                    if (atIdxC == atIdxB) continue;
                                    if (atIdxC == atIdx) continue;
                                    if (ConnMatrix[atIdxC][atIdxB] != 0) {

                                        if ( (!CurMol.getAtom(atIdxC).getFlag(CDKConstants.ISAROMATIC)) &&
                                                (!curRing.contains(CurMol.getAtom(atIdxC))) &&
                                                (ConnMatrix[atIdxC][atIdxC] == 6) ){
                                            int curVD = 0;
                                            int Odbl = 0;
                                            for (int k = 0; k < nAtoms; k++) {
                                                if (k == atIdxC) continue;
                                                if (ConnMatrix[atIdxC][k] != 0) {
                                                    curVD++;
                                                    if ( (ConnMatrix[atIdxC][k] == 2) && ConnMatrix[k][k] == 8)
                                                        Odbl++;
                                                }
                                            }
                                            if ( (curVD == 3) && (Odbl == 1) )
                                                FoundCO = atIdxC;
                                        }
                                    }
                                }

                                if ( (FoundOH > -1) && (FoundNH2_b > -1) )
                                    return true;
                                if ( (FoundOH > -1) && (FoundCO > -1) )
                                    return true;
                                if ( (FoundNH2 > -1) && (FoundCO > -1) )
                                    return true;
                            }
                        }
                    }
                }
            }

            // search NH2 - OH on different fused rings, and quinoline-type OH/NH2 -N=
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));

                int FoundNH2 = Find_Attached_NH2(atIdx);
                int FoundOH = Find_Attached_OH(atIdx);

                if ( (FoundNH2 > -1) || (FoundOH > -1)) {

                    // cycles on connected atoms to find one belonging to 2 different rings
                    for (int atIdxB = 0; atIdxB < nAtoms; atIdxB++) {
                        if (atIdx == atIdxB) continue;
                        if (ConnMatrix[atIdx][atIdxB] != 0) {

                            int countRings = 0;
                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                if (MolRings.getAtomContainer(rr).contains(CurMol.getAtom(atIdxB)))
                                    countRings++;
                            }
                            if (countRings != 2)
                                continue;

                            for (int atIdxC = 0; atIdxC < nAtoms; atIdxC++) {
                                if (atIdxC == atIdxB) continue;
                                if (atIdxC == atIdx) continue;
                                if ( (ConnMatrix[atIdxC][atIdxB] != 0) && (!curRing.contains(CurMol.getAtom(atIdxC)))) {
                                    // connected atoms at topo distance = 2 from atIdx, in a different but fused ring

                                    int FoundNH2_b = Find_Attached_NH2(atIdxC);
                                    int FoundOH_b = Find_Attached_OH(atIdxC);

                                    if ( (FoundNH2 > -1) && (FoundOH_b > -1) )
                                        return true;
                                    if ( (FoundOH > -1) && (FoundNH2_b > -1) )
                                        return true;

                                    // check quinolines
                                    if (ConnMatrix[atIdxC][atIdxC] == 7) {
                                        int nVD = 0, nSingle = 0, nDbl = 0;
                                        for (int k = 0; k < nAtoms; k++) {
                                            if (k == atIdxC) continue;
                                            if (ConnMatrix[atIdxC][k] != 0) {
                                                nVD++;
                                                if (ConnMatrix[atIdxC][k] == 1) nSingle++;
                                                if (ConnMatrix[atIdxC][k] == 2) nDbl++;
                                            }
                                        }
                                        if ( (nVD==2) && (nSingle==1) && (nDbl==2)) {
                                            // check ring sizes
                                            int RingSize1 = curRing.getAtomCount();
                                            int RingSize2 = 0;
                                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                                if (MolRings.getAtomContainer(rr).contains(CurMol.getAtom(atIdxC))) {
                                                    RingSize2 = MolRings.getAtomContainer(rr).getAtomCount();
                                                    break;
                                                }
                                            }
                                            if ( (RingSize1 == 6) && (RingSize2 == 6) )
                                                return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return false;
    }


    // Return atom index of attached NH2, or -1 if not found
    private int Find_Attached_NH2(int atIdx) {

        for (int i = 0; i < nAtoms; i++) {
            if (atIdx == i) continue;
            if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) continue;
            if ((ConnMatrix[atIdx][i] != 0)) {
                if (ConnMatrix[i][i] == 7) {
                    int nH;
                    try {
                        nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                    } catch (Exception e) {
                        nH = 0;
                    }
                    if (nH == 2)
                        return i;
                }
            }
        }

        return -1;
    }


    // Return atom index of attached OH, or -1 if not found
    private int Find_Attached_OH(int atIdx) {

        for (int i = 0; i < nAtoms; i++) {
            if (atIdx == i) continue;
            if (CurMol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) continue;
            if ((ConnMatrix[atIdx][i] != 0)) {
                if (ConnMatrix[i][i] == 8) {
                    int nH;
                    try {
                        nH = CurMol.getAtom(i).getImplicitHydrogenCount();
                    } catch (Exception e) {
                        nH = 0;
                    }
                    if (nH == 1)
                        return i;
                }
            }
        }

        return -1;
    }


    private double Check_NCS(IAtomContainer mol, int atom) {

        // Starting from S atom

        for (int i=0; i<nAtoms; i++) {
            if (i==atom) continue;
            if ( (ConnMatrix[i][atom]>0) && (AreEqual(GetName(mol,i),"C")) ) {
                for (int j=0; j<nAtoms; j++) {
                    if (j==i) continue;
                    if (AreEqual(GetName(mol,j),"N")) {
                        if (ConnMatrix[i][j]==2)
                            return 1;       // -N=C=S
                        if (ConnMatrix[i][j]==3)
                            return 0.5;     // N#C-S-
                    }
                }
            }
        }

        return 0;
    }


    private double Check_NCS() {

        double NCS = 0;

        for (int i=0; i<nAtoms; i++) {

            // Starting from C atom
            if (ConnMatrix[i][i] != 6) continue;

            int VD = 0;
            double bondS = 0;
            double bondN = 0;
            for (int j=0; j<nAtoms; j++) {
                if (i==j) continue;
                if (ConnMatrix[i][j] != 0) {
                    VD++;
                    if (ConnMatrix[i][j] == 7)
                        bondN = ConnMatrix[i][j];
                    if (ConnMatrix[i][j] == 16)
                        bondS = ConnMatrix[i][j];
                }
            }

            if (VD==2) {
                if ( (bondS == 2) && (bondN == 2) )
                    NCS += 1;
                if ( (bondS == 1) && (bondN == 3) )
                    NCS += 0.5;
            }

        }

        return NCS;
    }



    private int Check_BLM(IAtomContainer mol, int atom) {

        // Starting from N atom

        if (MolRings.contains(mol.getAtom(atom))) {
//            IRing R = (IRing)MolRings.getRings(mol.getAtom(atom)).getAtomContainer(0);
            for (IAtomContainer ac : MolRings.getRings(mol.getAtom(atom)).atomContainers() )  {
                IRing R = (IRing)ac;
                if (R.getAtomCount()==4) {
                    int N=0, C=0, CdblO=0, Other=0;
                    for (int i=0; i<4; i++) {
                        if (AreEqual(R.getAtom(i).getSymbol(),"N"))
                            N++;
                        else if (AreEqual(R.getAtom(i).getSymbol(),"C")) {
                            C++;
                            for (int k=0; k<nAtoms; k++)
                                if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[mol.indexOf(R.getAtom(i))][k]==2) ) {
                                    CdblO++;
                                }
                        } else
                            Other++;
                    }
                    if ( (N==1) && (C==3) && (CdblO==1) )
                        return 1;
                }
            }

        }

        return 0;
    }



    private int Calculate_POL() {

        int nPolSub = 0;

        for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
            IRing curRing = (IRing) MolRings.getAtomContainer(r);

            if (!IsRingPSA(curRing))
                continue;

            // check all ring atoms for substituent
            for (int i=0; i<curRing.getAtomCount(); i++) {
                int atIdx = CurMol.indexOf(curRing.getAtom(i));
                if (ConnMatrix[atIdx][atIdx] == 6) {

                    for (int j=0; j<nAtoms; j++) {
                        if (atIdx == j) continue;
                        if ( (ConnMatrix[atIdx][j] != 0) ) {

                            // check if in same ring or fused ring
                            boolean sameRing = false;
                            for (int rr=0; rr<MolRings.getAtomContainerCount(); rr++) {
                                IRing bufRing = (IRing) MolRings.getAtomContainer(rr);
                                // Dragon seems not to consider fused rings if they are not PSA
                                if (!IsRingPSA(bufRing))
                                    continue;
                                if  ( (bufRing.contains(CurMol.getAtom(j))) && (bufRing.contains(CurMol.getAtom(atIdx))) ) {
                                    sameRing = true;
                                    break;
                                }
                            }
                            if (sameRing)
                                continue;

                            if (ConnMatrix[j][j] != 6) {
                                // all hetero atoms considered pol sub
                                nPolSub++;
                            } else {
                                int nHet = 0;
                                int cVD = 0;
                                for (int k=0; k<nAtoms; k++) {
                                    if (k==j) continue;
                                    if (ConnMatrix[k][j] != 0) {
                                        cVD++;
                                        if (ConnMatrix[k][k] != 6)
                                            nHet++;
                                    }
                                }

                                int nH;
                                try {
                                    nH = CurMol.getAtom(j).getImplicitHydrogenCount();
                                } catch (Exception e) {
                                    nH = 0;
                                }
                                cVD += nH;

                                if ( (cVD == 2) && (nHet >= 1) )
                                    nPolSub++;
                                if ( (cVD == 3) && (nHet >= 1) )
                                    nPolSub++;
                                if ( (cVD == 4) && (nHet >= 2) )
                                    nPolSub++;
                            }
                        }
                    }
                }
            }

        }

        return nPolSub;
    }


    private boolean IsRingPSA(IRing curRing) {

        // consider aromatic PSA rings as in Dragon
        // Dragon seems NOT to consider if the ring structure is aromatic
        boolean AromaticPSA = true;
        for (int i=0; i<curRing.getAtomCount(); i++) {
            int atIdx = CurMol.indexOf(curRing.getAtom(i));
//                if (!CurMol.getAtom(atIdx).getFlag(CDKConstants.ISAROMATIC)) {
//                    AromaticPSA = false; break;
//                }
            int VD = 0;
            for (int j=0; j<nAtoms; j++) {
                if (atIdx == j) continue;
                if (ConnMatrix[atIdx][j] != 0)
                    VD++;
            }
            int nH;
            try {
                nH = CurMol.getAtom(atIdx).getImplicitHydrogenCount();
            } catch (Exception e) {
                nH = 0;
            }

            VD+= nH;
            if (VD > 3) {
                AromaticPSA = false;
                break;
            }
        }
        return AromaticPSA;
    }


    private int Calculate_RNG(IAtomContainer mol) {

        for (int i=0; i<MolRings.getAtomContainerCount(); i++) {

            IRing R = (IRing)MolRings.getAtomContainer(i);
            boolean SkipRing = false;

            // Skips benzene
            if (R.getAtomCount()==6)
                if (MoleculeUtilities.IsRingAromatic(R)) {
                    SkipRing = true;
                    for (int j=0; j<R.getAtomCount(); j++)
                        if (!(AreEqual(R.getAtom(j).getSymbol(),"C")))
                            SkipRing = false;
                }

            //
            if (!IsRingPSA(R))
                return 1;

            // Skips rings fused with benzene
            if (!SkipRing) {
                for (int at=0; at<R.getAtomCount(); at++) {
                    IAtom At = R.getAtom(at);
                    for (int r=0; r<MolRings.getAtomContainerCount(); r++) {
                        if (r == i) continue;
                        IRing curR = (IRing)MolRings.getAtomContainer(r);
                        if (curR.contains(At)) {
                            // Checks if this ring is benzene
                            if ( (curR.getAtomCount() ==6 ) &&
                                    (MoleculeUtilities.IsRingAromatic(curR)) ) {
                                SkipRing = true;
                                break;
                            }
                        }
                    }
                    if (SkipRing) break;
                }
            }

            if (!SkipRing)
                return 1;
        }

        return 0;
    }



    private int Calculate_ALK(IAtomContainer mol) {

        boolean OnlyC;
        int nC=0;

        // Check if it is a hydrocarbon compound
        OnlyC = true;
        for (int i=0; i<nAtoms; i++) {
            if  (ConnMatrix[i][i] == 6) nC++;
            else OnlyC = false;
        }

        // Hydrocarbon compound
        if (OnlyC) {
            boolean HasAromaticBonds = false;
            for (int i=0; i<nAtoms; i++)
                for (int j=0; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j] == 1.5) {
                        HasAromaticBonds = true;
                        break;
                    }
                }

            if (HasAromaticBonds)
                return 0;

            // Checks alkane and alkene
            int nSingle=0, nDouble=0, nBonds=0;
            for (int i=0; i<(nAtoms-1); i++)
                for (int j=i; j<nAtoms; j++) {
                    if (i==j) continue;
                    if (ConnMatrix[i][j]>0) {
                        nBonds++;
                        if (ConnMatrix[i][j] == 1) nSingle++;
                        if (ConnMatrix[i][j] == 2) nDouble++;
                    }
                }
            if (nBonds == nSingle)
                return 1; // alkane

        }

        // Checks for hydrocarbon chains with at least 7 C
        if (nC>6) {

            for (int i=0; i<nAtoms; i++)
                if  (AreEqual(GetName(mol,i),"C")) {

                    IAtom chainAt = mol.getAtom(i);
                    int nH=0;
                    if (chainAt.getImplicitHydrogenCount()!=null)
                        nH = chainAt.getImplicitHydrogenCount();
                    if (nH!=3)
                        continue; // shall start from terminal CH3 atom

                    for (int j=0; j<nAtoms; j++)
                        if ((AreEqual(GetName(mol,j),"C"))&& (TopoMatrix[i][j]>6)) {

                            boolean ChainFound=true;

                            ShortestPaths shortestPaths = new ShortestPaths(mol, mol.getAtom(i));
                            List<IAtom> Path = Arrays.asList(shortestPaths.atomsTo(mol.getAtom(j)));

                            // DEPRECATED METHOD
//                            List<IAtom> Path = PathTools.getShortestPath(mol, mol.getAtom(i), mol.getAtom(j));

                            for (int k=1; k<Path.size(); k++) {
                                IAtom CurAtom = Path.get(k);
                                nH=0;
                                if (CurAtom.getImplicitHydrogenCount()!=null)
                                    nH = CurAtom.getImplicitHydrogenCount();
                                if ( (!(CurAtom.getSymbol().equalsIgnoreCase("C"))) ||
                                        (MolRings.contains(CurAtom)) ||
                                        (nH!=2) ) {
                                    // If atom is no Carbon or is not a chain, breaks
                                    ChainFound = false;
                                    break;
                                }
                            }
                            if (ChainFound)
                                return 1;
                        }
                }
        }


        return 0;
    }



    private int Calculate_AMP(IAtomContainer mol) {

        int count = 0;

        // Checks for COOH
        for (int i=0; i<nAtoms; i++) {

            boolean isC_COOH = false;
            if (!(AreEqual(GetName(mol,i),"C")))
                for (int j=0; j<nAtoms; j++) {
                    if ( (AreEqual(GetName(mol,j),"C")) && (ConnMatrix[i][j]==1) ) {
                        int Osingle=0, Odbl=0;
                        for (int k=0; k<nAtoms; k++) {
                            if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[j][k]==1) )
                                Osingle++;
                            if ( (AreEqual(GetName(mol,k),"O")) && (ConnMatrix[j][k]==2) )
                                Odbl++;
                        }
                        if ( (Osingle==1) && (Odbl==1) )
                            isC_COOH = true;
                    }
                }

            if (isC_COOH) {

                Boolean Pyridine=false, Amine=false;

                if (MolRings.contains(mol.getAtom(i))) {

                    // C atom is in ring, checks for Pyridine or aminobenzoic

                    boolean InvalidStruct = false;
                    IRing R = (IRing)MolRings.getRings(mol.getAtom(i)).getAtomContainer(0);
                    if ( (R.getAtomCount()==6) && (MoleculeUtilities.IsRingAromatic(R)) ) {
                        for (int j=0; j<6; j++) {
                            int Branches = 0;
                            int atomNum = mol.indexOf(R.getAtom(j));
                            if (atomNum!=i) {
                                if (AreEqual(GetName(mol,j),"C")) {
                                    for (int k=0; k<nAtoms; k++) {
                                        if ( (ConnMatrix[k][atomNum]>0) && (!(R.contains(mol.getAtom(k)))) ) {
                                            if ( (AreEqual(GetName(mol,k),"N")) && (GetHydrogens(mol, k)==2)) {
                                                Branches++;
                                                Amine = true;
                                            }
                                        }
                                    }
                                    if ((Branches>0)&&(!Amine)) {
                                        InvalidStruct = true;
                                        break;
                                    }
                                }
                                if (AreEqual(GetName(mol,j),"N")) {
                                    for (int k=0; k<nAtoms; k++) {
                                        if ( (ConnMatrix[k][atomNum]>0) && (!(R.contains(mol.getAtom(k)))) ) {
                                            InvalidStruct = true;
                                            break;
                                        }
                                    }
                                    Pyridine = true;
                                }
                            }
                        }

                        if ( (!InvalidStruct) && (Pyridine||Amine) )
                            count += 0.5;
                    }

                } else {

                    // No cyclic structure, checks for AA

                    int NH2=0, C=0, R=0;
                    for (int k=0; k<nAtoms; k++)
                        if (ConnMatrix[k][i]>0) {
                            if (AreEqual(GetName(mol,k),"N"))
                                if (GetHydrogens(mol, k)==2)
                                    NH2++;
                                else
                                    R++;
                            else if (AreEqual(GetName(mol,k),"C"))
                                C++;
                            else
                                R++;
                        }
                    if (NH2==1)
                        if ( ((C==1)&&(R==1)) || (C==2) )
                            count += 1;

                }

            }

        }

        return count;
    }



    ////// Utilities for this class

    private String GetName(IAtomContainer mol, int atom) {
        return mol.getAtom(atom).getSymbol();
    }

    private int GetHydrogens(IAtomContainer mol, int atom) {

        int nH=0;
        if (mol.getAtom(atom).getImplicitHydrogenCount()!=null)
            nH = mol.getAtom(atom).getImplicitHydrogenCount();

        return nH;
    }

    private boolean AreEqual(String AtomA, String AtomB) {
        return (AtomA.equalsIgnoreCase(AtomB));
    }


    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        MLogP block = new MLogP();
        block.CloneDetailsFrom(this);
        return block;
    }
}
