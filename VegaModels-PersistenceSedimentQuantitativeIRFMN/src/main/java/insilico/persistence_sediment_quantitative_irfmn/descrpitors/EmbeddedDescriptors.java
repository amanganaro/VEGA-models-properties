package insilico.persistence_sediment_quantitative_irfmn.descrpitors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.persistence_sediment_quantitative_irfmn.descrpitors.weights.*;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;

@Log4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double MlogP = MISSING_VALUE;
    private double B3C_C = MISSING_VALUE;
    private double C_002 = MISSING_VALUE;
    private double CATS2D_3_DL = MISSING_VALUE;
    private double B4O_Cl = MISSING_VALUE;
    private double Gmax = MISSING_VALUE;
    private double EEig10dm = MISSING_VALUE;
    private double BEL8p = MISSING_VALUE;

    private double[][] ConnAugMatrixCats2D;
    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";
    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };


    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateMLogP(Mol);
        CalculateACF(Mol);
        CalculateTopoDistances(Mol);
        CalculateCATS2D(Mol);
        CalculateGmax(Mol);
        CalculateEEig(Mol);
        CalculateBEL(Mol);
    }

    private void CalculateBEL(InsilicoMolecule Mol){
        this.setBEL8p(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setBEL8p(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            setBEL8p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        double[] w = Polarizability.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            if(i == 7){
                double  valL;
                if (i>(eigenvalues.length-1)) {
                    valL = 0;
                } else {
                    valL = eigenvalues[i];
                }
                this.setBEL8p(valL);
                return;
            }
        }
    }

    private void CalculateEEig(InsilicoMolecule Mol){
        this.setEEig10dm(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setEEig10dm(MISSING_VALUE);
            return;
        }

        if (m.getAtomCount() < 2) {
            this.setEEig10dm(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setEEig10dm(MISSING_VALUE);
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeDipoleMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++)
            for (int j=0; j<EdgeAdjMat[0].length; j++)
                EdgeDipoleMat[i][j] = EdgeAdjMat[i][j][0];

        for (int i=0; i<m.getBondCount(); i++) {
            for (int j=0; j<m.getBondCount(); j++) {
                if (EdgeDipoleMat[i][j] != 0) {
                    IAtom a = m.getBond(i).getAtom(0);
                    IAtom b = m.getBond(i).getAtom(1);
                    double CurVal = GetDipoleMoment(m, a, b);
                    if (CurVal == 0)
                        CurVal = GetDipoleMoment(m, b, a);
                    EdgeDipoleMat[i][j] = CurVal;
                }
            }
        }

        DataMatrix = new Matrix(EdgeDipoleMat);

        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);


        // EEig
        for (int i=1; i<=15; i++) {
            if(i == 10) {
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    this.setEEig10dm(eigenvalues[idx]);
                else
                    this.setEEig10dm(0);
                return;
            }
        }
    }

    private void CalculateGmax(InsilicoMolecule Mol){
        this.setGmax(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setGmax(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            this.setGmax(MISSING_VALUE);
            return;
        }

        double Gmax= MISSING_VALUE;

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Maximum and minimum Estate/HEstate
            Gmax = (Gmax == Descriptor.MISSING_VALUE) ? es.getEState()[at] : (Math.max(es.getEState()[at], Gmax));
        }

        this.setGmax(Gmax);

    }

    private void CalculateCATS2D(InsilicoMolecule Mol){
        this.setCATS2D_3_DL(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }
        try {
            ConnAugMatrixCats2D = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_3_DL(MISSING_VALUE);
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (int d=0; d<AtomCouples.length; d++) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i=0; i<nSK; i++) {
                if ( isIn(CatsTypes[i], AtomCouples[d][0]) ) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if ( isIn(CatsTypes[j], AtomCouples[d][1]) ) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i=0; i<desc.length; i++)
                if(i == 3)
                    this.setCATS2D_3_DL(desc[i]);

        }


    }

    private void CalculateTopoDistances(InsilicoMolecule Mol){
        this.setB3C_C(0); this.setB4O_Cl(0);

        String[][] AtomCouplesCC = {
            {"C", "C"},
        };

        String[][] AtomCouplesOCl = {
                {"O", "Cl"},
        };


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setB3C_C(MISSING_VALUE); this.setB4O_Cl(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setB3C_C(MISSING_VALUE); this.setB4O_Cl(MISSING_VALUE);
            return;
        }

        int[] descB = new int[10];
        Arrays.fill(descB, 0);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesCC[0][0])) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesCC[0][1])) {

                        // B (presence of pair) and F (number of couples)
                        if (TopoMat[i][j] <= 10) {
                            descB[TopoMat[i][j]-1] = 1;
                        }

                    }
                }
            }
        }


        for (int i=0; i<descB.length; i++) {
            if(i == 2)
                this.setB3C_C(descB[i]);
        }

        descB = new int[10];
        Arrays.fill(descB, 0);

        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesOCl[0][0])) {
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesOCl[0][1])) {

                        // B (presence of pair) and F (number of couples)
                        if (TopoMat[i][j] <= 10) {
                            descB[TopoMat[i][j]-1] = 1;
                        }

                    }
                }
            }
        }


        for (int i=0; i<descB.length; i++) {
            if(i == 3)
                this.setB4O_Cl(descB[i]);
        }





    }

    private void CalculateACF(InsilicoMolecule Mol){
        this.setC_002(0);
        DescriptorBlock descriptorBlock = new ACF();
        try {
            descriptorBlock.Calculate(Mol);
            this.setC_002(descriptorBlock.GetByName("C-002").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setC_002(MISSING_VALUE);
        }
    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        this.setMlogP(0);
        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        this.setMlogP(descriptorMLogP.getMLogP());
    }

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }

    /**
     * Sets CATS 2D atom types for each atom, as a list of string containing
     * all matching types for each atom.
     **/
    private ArrayList<String>[] setCatsAtomType(IAtomContainer m) {

        int nSK = m.getAtomCount();
        ArrayList<String>[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false;

            // Definition of CATS types
            //
            // A: O, N without H
            // N: [+], NH2
            // P: [-], COOH, POOH, SOOH

            // Hydrogens
            int H = 0;
            try {
                H = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) { }

            // [+]
            if (CurAt.getFormalCharge() > 0) {

                boolean NpOm = false;
                if (ConnAugMatrixCats2D[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrixCats2D[i][j]==1) {
                            if (ConnAugMatrixCats2D[j][j] == 8) {
                                IAtom Oxy = m.getAtom(j);
                                if (Oxy.getFormalCharge()!=0)
                                    NpOm = true;
                            }
                        }
                    }
                }

                if (!NpOm)
                    tP = true;
            }

            // [-]
            if (CurAt.getFormalCharge() < 0)
                tN = true;

            // O
            if (CurAt.getSymbol().equalsIgnoreCase("O")) {
                tA = true;

                if ( (CurAt.getFormalCharge() == 0) && (H == 1))
                    tD = true;

            }

            // N (NH2 and N without H)
            if (CurAt.getSymbol().equalsIgnoreCase("N")) {

                int nSglBnd = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if (ConnAugMatrixCats2D[i][j] == 1)
                            nSglBnd++;
                        else
                            nOtherBnd++;
                    }
                }

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (CurAt.getFormalCharge() == 0) &&( (H == 1) || (H ==2) ) )
                    tD = true;

            }

            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                    (CurAt.getFormalCharge() == 0) )  {

                int nSglBnd = 0, nDblO = 0, nSglOH = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if (ConnAugMatrixCats2D[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrixCats2D[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrixCats2D[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrixCats2D[i][j] == 2) && (ConnAugMatrixCats2D[j][j] == 8) )
                                nDblO++;
                            else
                                nOtherBnd++;
                        }
                    }
                }

                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherBnd == 0) )
                    tN = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                boolean connOnlyToSingleC = true;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if ( (ConnAugMatrixCats2D[j][j] != 6) || (ConnAugMatrixCats2D[i][j] > 1.5) ) {
                            connOnlyToSingleC = false;
                            break;
                        }
                    }
                }
                if (connOnlyToSingleC)
                    tL = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                boolean connOnlyToSingleC = true;
                int nSingleC = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCats2D[i][j]>0) {
                        if ( (ConnAugMatrixCats2D[j][j] != 6) || (ConnAugMatrixCats2D[i][j] != 1) ) {
                            connOnlyToSingleC = false;
                            break;
                        } else {
                            nSingleC++;
                        }
                    }
                }
                if ( (connOnlyToSingleC) && (nSingleC == 2) )
                    tL = true;
            }


            // Sets final types
            if (tD) AtomTypes[i].add(TYPE_D);
            if (tL) AtomTypes[i].add(TYPE_L);

        }

        return AtomTypes;
    }

    private double GetDipoleMoment(IAtomContainer CurMol, IAtom at1, IAtom at2) {

        String a = at1.getSymbol();
        String b = at2.getSymbol();

        // C - something
        if (a.equalsIgnoreCase("C")) {

            // C-F
            if (b.equalsIgnoreCase("F")) {
                return 1.51;
            }

            // C-Cl , C(Cl)-Cl , C(Cl)(Cl)-Cl
            if (b.equalsIgnoreCase("Cl")) {
                int nCl=0;
                for (IAtom at : CurMol.getConnectedAtomsList(at1)) {
                    if (at.getSymbol().equalsIgnoreCase("Cl"))
                        nCl++;
                }
                if (nCl==1)
                    return 1.56;
                if (nCl==2)
                    return 1.20;
                if (nCl==3)
                    return 0.83;
            }

            // C-Br
            if (b.equalsIgnoreCase("Br")) {
                return 1.48;
            }

            // C-I
            if (b.equalsIgnoreCase("I")) {
                return 1.29;
            }

            // C-N , C=N , C#N
            if (b.equalsIgnoreCase("N")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.4;
                if (ord == IBond.Order.DOUBLE)
                    return 0.9;
                if (ord == IBond.Order.TRIPLE)
                    return 3.6;
            }

            // C-O , C=O
            if (b.equalsIgnoreCase("O")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.86;
                if (ord == IBond.Order.DOUBLE)
                    return 2.4;
            }

            // C-S , C=S
            if (b.equalsIgnoreCase("S")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 2.95;
                if (ord == IBond.Order.DOUBLE)
                    return 2.8;
            }

        }


        // N-O , N-[O-] , N=O
        if ((a.equalsIgnoreCase("N")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nH=0;
            try { nH = at2.getImplicitHydrogenCount(); } catch (Exception e) {}
            int nConn = CurMol.getConnectedBondsCount(at2) + nH;
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 3.2;
            if ((ord == IBond.Order.SINGLE) && (nConn==1))
                return 2.0;
//                return 0.3;
            if ((ord == IBond.Order.DOUBLE) && (nConn==1))
                return 2.0;
        }


        // S-[O-]
        if ((a.equalsIgnoreCase("S")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nConn = CurMol.getConnectedBondsCount(at2);
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 2.9;
        }


        // C(*)(*)-C(*)(*)(*) , C(*)(*)-C , CC(*)(*)(*)
        if ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("C"))) {

            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;

            int nH1=0, nH2=0;
            try {
                nH1 = at1.getImplicitHydrogenCount();
            } catch (Exception E) {}
            try {
                nH2 = at2.getImplicitHydrogenCount();
            } catch (Exception E) {}

            int nConn1 = CurMol.getConnectedBondsCount(at1) + nH1;
            int nConn2 = CurMol.getConnectedBondsCount(at2) + nH2;

            if ((nConn1==3) && (nConn2==4))
                return 0.68;
            if ((nConn1==3) && (nConn2==2))
                return 1.15;
            if ((nConn1==2) && (nConn2==4))
                return 1.48;
        }

        return 0;
    }


}
