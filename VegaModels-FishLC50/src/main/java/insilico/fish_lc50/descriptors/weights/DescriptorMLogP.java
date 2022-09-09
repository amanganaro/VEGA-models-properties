package insilico.fish_lc50.descriptors.weights;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.log4j.Log4j;
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


public class DescriptorMLogP {
    private static final Logger log = LogManager.getLogger(DescriptorMLogP.class);

    private double MISSING_VALUE = -999;
    public double MLogP = MISSING_VALUE;
    private int nAtoms;
    private double[][] ConnMatrix;
    private int[][] TopoMatrix;
    private RingSet MolRings;

    public DescriptorMLogP(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    public void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateMLogP(Mol);
    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        MLogP = 0;
        IAtomContainer CurMol;

        try {
            CurMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MLogP = MISSING_VALUE;
            return;
        }

        // Gets matrices
        try {
            ConnMatrix = Mol.GetMatrixConnectionAugmented();
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MLogP = MISSING_VALUE;
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
//        SSSRFinder rf = new SSSRFinder(CurMol);
//        MolRings = (RingSet) rf.findSSSR();
        try {
            MolRings = Mol.GetSSSR();
        } catch (InvalidMoleculeException e) {
            MLogP = MISSING_VALUE;
            return;
        }

        POL = Calculate_POL(CurMol);
        RNG = Calculate_RNG(CurMol);
        ALK = Calculate_ALK(CurMol);
        AMP = Calculate_AMP(CurMol);


        // UB
        //
        // TODO - check problems with bond order with aromatic bonds
        //
        int nArBonds = 0;
        for (int i=0; i<CurMol.getBondCount(); i++) {
            if (CurMol.getBond(i).getFlag(CDKConstants.ISAROMATIC)) {
                nArBonds++;
                continue;
            }
            IBond.Order bo = CurMol.getBond(i).getOrder();
            if ((bo == IBond.Order.DOUBLE) || (bo == IBond.Order.TRIPLE))
                UB++;
        }
        // approximation of double bonds in aromatic rings
        UB += Math.floor((double)nArBonds / 2.00);


        for (int i=0; i<nAtoms; i++) {

            // Checks on single atom type

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("C")) {

                CX++;

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("O")) {

                NO++;
                PRX += Check_PRX(CurMol, i);
                if (HB == 0)
                    HB = Check_HB(CurMol, i);

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("N")) {

                NO++;
                PRX += Check_PRX(CurMol, i);

                // Checks for NO2 for UB and >N< for QN and NO2
                int nOdbl=0;
                int nOminus=0;
                int nVD=0;
                for (int j=0; j<nAtoms; j++) {
                    if (j==i) continue;
                    if (ConnMatrix[i][j]>0) {
                        nVD++;
                        if ((ConnMatrix[i][j]==2) && (ConnMatrix[j][j]==8))
                            nOdbl++;
                        if ((ConnMatrix[i][j]==1) && (ConnMatrix[j][j]==8)) {
                            int OVD=0;
                            for (int k=0; k<nAtoms; k++) {
                                if (k==j) continue;
                                if (ConnMatrix[j][k]>0) OVD++;
                            }
                            if (OVD==1)
                                nOminus++;
                        }
                    }
                }

                // Quaternary nitrogen or N-oxide
                if (nVD==4) {
                    if (nOminus>0)
                        QN += 0.5;
                    else
                        QN += 1;
                }

                // Nitro group
                if ((nVD==3) && (nOdbl==1) && (nOminus==1)) {
                    // note: due to preprocessing, nitro groups are always
                    // set in O=[N+][O-] form
                    NO2++;
                }


                if (HB == 0)
                    HB = Check_HB(CurMol, i);

                if (BLM == 0)
                    BLM = Check_BLM(CurMol, i);

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("S")) {

                NCS += Check_NCS(CurMol, i);

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("F")) {

                CX += 0.5;

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("Cl")) {

                CX++;

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("Br")) {

                CX += 1.5;

                continue;
            }

            if (CurMol.getAtom(i).getSymbol().equalsIgnoreCase("I")) {

                CX += 2;
            }

        }

        PRX = PRX/2;

        // Calculation of LogP
        MLogP = -1.041 + 1.244 * Math.pow(CX,0.6) - 1.017 * Math.pow(NO, 0.9) +
                0.406 * PRX + -0.145 * Math.pow(UB, 0.8) + 0.511 * HB +
                0.268 * POL - 2.215 * AMP + 0.912 * ALK -0.392 * RNG -
                3.684 * QN + 0.474 * NO2 + 1.582 * NCS + 0.773 * BLM;

        // Sets final descriptor value
    }




    private int Check_PRX(IAtomContainer mol, int atom) {

        int count=0;

        if (!( (AreEqual(GetName(mol, atom),"O")) || (AreEqual(GetName(mol, atom),"N")) ))
            return 0;

        for (int i=0; i<nAtoms; i++) {

            if ( (i==atom) )
                continue;

            // X-Y
            if ((TopoMatrix[i][atom]==1) && ((ConnMatrix[i][i]==8)||(ConnMatrix[i][i]==7))) {
                count+=2;



            }

            // X-A-Y
            if ((TopoMatrix[i][atom]==2)  && ((ConnMatrix[i][i]==8)||(ConnMatrix[i][i]==7))) {

                ShortestPaths shortestPaths = new ShortestPaths(mol, mol.getAtom(atom));
                List<IAtom> Path = Arrays.asList(shortestPaths.atomsTo(mol.getAtom(i)));

                // DEPRECATED METHOD
//                List<IAtom> Path = PathTools.getShortestPath(mol, mol.getAtom(atom), mol.getAtom(i));

                String A_Atom = Path.get(1).getSymbol();
                if ((AreEqual(A_Atom,"C"))||(AreEqual(A_Atom,"S"))||(AreEqual(A_Atom,"P"))) {
                    count+=1;


                    int MidAtom = mol.indexOf(Path.get(1));
                    int midVD = 0;
                    for (int k=0; k<nAtoms; k++) {
                        if (k==MidAtom) continue;
                        if (ConnMatrix[MidAtom][k]>0) midVD++;
                    }

                    if (midVD==2) {

                        int atomVD = 0;
                        for (int k=0; k<nAtoms; k++) {
                            if (k==atom) continue;
                            if (ConnMatrix[atom][k]>0) atomVD++;
                        }
                        int iVD = 0;
                        for (int k=0; k<nAtoms; k++) {
                            if (k==i) continue;
                            if (ConnMatrix[i][k]>0) iVD++;
                        }

                        int midOsingle=0, midOdouble=0, midNsingle=0;
                        // atom (X)
                        if (ConnMatrix[MidAtom][atom]==1) {
                            if ((ConnMatrix[atom][atom]==7)&&(atomVD==3)) midNsingle++;
                            if ((ConnMatrix[atom][atom]==8)&&(atomVD==2)) midOsingle++;
                        }
                        if (ConnMatrix[MidAtom][atom]==2) {
                            if ((ConnMatrix[atom][atom]==8)&&(atomVD==1)) midOdouble++;
                        }
                        // i (Y)
                        if (ConnMatrix[MidAtom][i]==1) {
                            if ((ConnMatrix[i][i]==7)&&(atomVD==3)) midNsingle++;
//                            if ((ConnMatrix[i][i]==8)&&(atomVD==2)) midOsingle++;
                            if (ConnMatrix[i][i]==8) midOsingle++;
                        }
                        if (ConnMatrix[MidAtom][i]==2) {
                            if ((ConnMatrix[i][i]==8)&&(atomVD==1)) midOdouble++;
                        }

                        boolean XaYvalue1=false;
                        if (midOsingle==2) XaYvalue1=true;
                        else if (midNsingle==2) XaYvalue1=true;
                        else if ((midOsingle==1)&&(midNsingle==1)) XaYvalue1=true;
                        else if (midOdouble == 1 && midOsingle == 1) XaYvalue1=true;

                        if (!(XaYvalue1))
                            count++;

                    } else if (midVD==3) {

                        // to be consistent with Dragon 5.5 also
                        // groups -R(=O)O have value 2.

                        int midOsingle=0, midOdouble=0;
                        // atom (X)
                        if (ConnMatrix[atom][atom]==8) {
                            if (ConnMatrix[MidAtom][atom]==1) midOsingle++;
                            if (ConnMatrix[MidAtom][atom]==2) midOdouble++;
                        }
                        // i (Y)
                        if (ConnMatrix[i][i]==8) {
                            if (ConnMatrix[MidAtom][i]==1) midOsingle++;
                            if (ConnMatrix[MidAtom][i]==2) midOdouble++;
                        }

                        if ((midOsingle==1)&&(midOdouble==1)) count++;

                    }

                }
            }
        }

        //
        // TODO
        // mancano le correction -CON< e -SO2N<
        //

        return count;
    }


    private int Check_HB(IAtomContainer mol, int atom) {

        boolean isOHorNH2 = false;

        for (int i=0; i<nAtoms; i++) {
            if (i==atom)
                continue;

            if (ConnMatrix[i][atom]>0) {
                if ( (ConnMatrix[i][atom]==1) && (!isOHorNH2) )
                    isOHorNH2 = true;
                else
                    isOHorNH2 = false;
            }
        }

        if (!isOHorNH2)
            return 0;

        for (int i=0; i<nAtoms; i++) {

            if (i==atom)
                continue;

            if ( (TopoMatrix[i][atom] == 1) && (mol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) ) {
                for (int j=0; j<nAtoms; j++) {
                    if ( (TopoMatrix[atom][j] == 3) && (TopoMatrix[atom][j] == 4) )
                        if ((AreEqual(GetName(mol,j),"O"))||(AreEqual(GetName(mol,j),"N")))
                            return 1;
                }
            }
        }

        return 0;
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



    private int Calculate_POL(IAtomContainer mol) {

        int count = 0;

        for (int i=0; i<MolRings.getAtomContainerCount(); i++) {
            IRing R = (IRing)MolRings.getAtomContainer(i);
//            if (AromaticityCalculator.isAromatic(R, mol))
            if (MoleculeUtilities.IsRingAromatic(R))
                for (int j=0; j<R.getAtomCount(); j++) {

                    int jj = mol.indexOf(R.getAtom(j));

                    for (int k=0; k<nAtoms; k++)
                        if ( (TopoMatrix[k][jj]==1) && (!(R.contains(mol.getAtom(k)))) &&
                                (!(MolRings.contains(mol.getAtom(k)))) ) {

                            boolean RootIsC=false, ChainLongerThan2=false;
                            int Other=0, sC=0, dC=0;

                            if (AreEqual(GetName(mol,k),"C"))
                                RootIsC = true;

                            for (int z=0; z<nAtoms; z++)
                                if (TopoMatrix[k][z]==1) {

                                    if (AreEqual(GetName(mol,z),"C")) {
                                        if (ConnMatrix[k][z]==1) sC++;
                                        else if (ConnMatrix[k][z]==2) dC++;
                                        else Other++;
                                    } else {
                                        Other++;
                                    }

                                    for (int zz=0; zz<nAtoms; zz++) {
                                        if (zz==z) continue;
                                        if (TopoMatrix[z][zz]==1)
                                            ChainLongerThan2 = true;
                                    }

                                }

                            boolean toSkip = false;

                            // skip excluded substituents
                            if (RootIsC) {
                                if ((Other==1)&&(dC==0)) toSkip = true;
                                if ((Other==0)&&(sC>0)&&(dC==0)) toSkip = true;
                                if ((Other==0)&&(dC==1)) toSkip = true;
                            }

                            if (!(toSkip)) {
                                count++;
                                if (count==4)
                                    return count;
                                break;
                            }

                        }

                }
        }

        return count;
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

            // Skips rings fised with benzene
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

}
