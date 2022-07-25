package insilico.persistence_quantative_water_irfmn.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.persistence_quantative_water_irfmn.descriptors.weights.*;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.ArrayList;
import java.util.Arrays;

@Log4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double MLogP = MISSING_VALUE;
    private double F1_C_C = MISSING_VALUE;
    private double CATS2D_3_DL = MISSING_VALUE;
    private double MATS8v = MISSING_VALUE;
    private double MATS1e = MISSING_VALUE;
    private double O_060 = MISSING_VALUE;
    private double MAXDP = MISSING_VALUE;
    private double C_008 = MISSING_VALUE;
    private double MATS6p = MISSING_VALUE;
    private double nR9 = MISSING_VALUE;
    private double F8_C_Cl = MISSING_VALUE;
    private double ATSC6e = MISSING_VALUE;
    private double GATS6v = MISSING_VALUE;
    private double F4_C_C = MISSING_VALUE;
    private double F3_N_Cl = MISSING_VALUE;
    private double MATS8m = MISSING_VALUE;

    private double[][] ConnAugMatrixCats2D;
    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";
    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };


    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateMLogP(Mol);
        CalculateACF(Mol);
        CalculateTopoDistances(Mol);
        CalculateCats2D(Mol);
        CalculateTopologicaEState(Mol);
        CalculateRings(Mol);
        CalculateAutocorrelation(Mol);
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol){
        this.setMATS1e(0); this.setMATS6p(0);
        this.setMATS8m(0); this.setMATS8v(0);
        this.setGATS6v(0); this.setATSC6e(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMATS1e(MISSING_VALUE); this.setMATS6p(MISSING_VALUE);
            this.setMATS8m(MISSING_VALUE); this.setMATS8v(MISSING_VALUE);
            this.setGATS6v(MISSING_VALUE); this.setATSC6e(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setMATS1e(MISSING_VALUE); this.setMATS6p(MISSING_VALUE);
            this.setMATS8m(MISSING_VALUE); this.setMATS8v(MISSING_VALUE);
            this.setGATS6v(MISSING_VALUE); this.setATSC6e(MISSING_VALUE);
            return;
        }
        
        int nSK = m.getAtomCount();

        // e
        double[] w = Electronegativity.getWeights(m);
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 1) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        // Sets descriptors
        this.setMATS1e(MoranAC);

        double  ACS=0;
        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 6) {
                    ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                }
        }

        ACS /= 2.0;
        this.setATSC6e(ACS);


        // p
        w = Polarizability.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 6) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        // Sets descriptors
        this.setMATS6p(MoranAC);

        //v
        w = VanDerWaals.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 8) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        // Sets descriptors
        this.setMATS8v(MoranAC);


        MoranAC=0;
        double GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 6) {
                    GearyAC += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                GearyAC = 0;
            } else {
                GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
            }
        }


        this.setGATS6v(GearyAC);


        // m
        w = Mass.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);
        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 8) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        // Sets descriptors
        this.setMATS8m(MoranAC);
    }

    private void CalculateRings(InsilicoMolecule Mol){
        this.setNR9(0);

        try {

            int nSizes = 9;
            int[] RingCount = new int[nSizes];
            int[] RingSize = new int[nSizes];
            for (int i=0; i<nSizes; i++) {
                RingSize[i] = 3 + i;
                RingCount[i] = 0;
            }

            IRingSet allRings = Mol.GetAllRings();
            for (IAtomContainer iAtomContainer : allRings.atomContainers()) {
                IRing ring = (IRing) iAtomContainer;
                for (int i = 0; i < nSizes; i++) {
                    if (ring.getAtomCount() == RingSize[i])
                        RingCount[i]++;
                }
            }

            for (int i=0; i<nSizes; i++)
                if(RingSize[i] == 9)
                    this.setNR9(RingCount[i]);

        } catch (Throwable e) {
            this.setNR9(MISSING_VALUE);
        }
    }

    private void CalculateTopologicaEState(InsilicoMolecule Mol){
        this.setMAXDP(0);

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMAXDP(MISSING_VALUE);
            return;
        }

        // Get weights
        EState ES;
        try {
            ES = new EState(Mol.GetStructure());
        } catch (Exception e) {
            log.warn(e.getMessage());
            this.setMAXDP(MISSING_VALUE);
            return;
        }
        double[] IStates = ES.getIS();

        // Get matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setMAXDP(MISSING_VALUE);
            return;
        }


        try {

            int nSK = curMol.getAtomCount();
            double maxDN = 0, maxDP = 0;

            for (int i=0; i<nSK; i++) {
                double Delta = 0;
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if ( (IStates[i]!= Descriptor.MISSING_VALUE) && (IStates[j]!= Descriptor.MISSING_VALUE) )
                        Delta += (IStates[i]-IStates[j]) / Math.pow(TopoMat[i][j]+1, 2);
                }
                if (Delta<0) {
                    if (Delta<maxDN)
                        maxDN = Delta;
                } else {
                    if (Delta>maxDP)
                        maxDP = Delta;
                }
            }


            this.setMAXDP(maxDP);

        } catch (Throwable e) {
            this.setMAXDP(MISSING_VALUE);
        }
    }

    private void CalculateCats2D(InsilicoMolecule Mol){
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
        this.setF1_C_C(0); this.setF3_N_Cl(0);
        this.setF4_C_C(0); this.setF8_C_Cl(0);

        String[][] AtomCouplesCC = {
                {"C", "C"},
        };

        String[][] AtomCouplesNCl = {
                {"N", "Cl"},
        };

        String[][] AtomCouplesCCl= {
                {"C", "Cl"},
        };

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setF1_C_C(MISSING_VALUE); this.setF3_N_Cl(MISSING_VALUE);
            this.setF4_C_C(MISSING_VALUE); this.setF8_C_Cl(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setF1_C_C(MISSING_VALUE); this.setF3_N_Cl(MISSING_VALUE);
            this.setF4_C_C(MISSING_VALUE); this.setF8_C_Cl(MISSING_VALUE);
            return;
        }

        //CC
        for (int d=0; d<AtomCouplesCC.length; d++) {

            int[] descF = new int[10];
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesCC[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesCC[d][1])) {

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= 10) {
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouplesCC[d][0].compareTo(AtomCouplesCC[d][1]) == 0) {
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i=0; i<descF.length; i++) {
                if(i == 0)
                    this.setF1_C_C(descF[i]);
                if(i == 3)
                    this.setF4_C_C(descF[i]);
            }
        }

        for (int d=0; d<AtomCouplesCC.length; d++) {

            int[] descF = new int[10];
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesCCl[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesCCl[d][1])) {

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= 10) {
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouplesCCl[d][0].compareTo(AtomCouplesCCl[d][1]) == 0) {
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i=0; i<descF.length; i++) {
                if(i == 7)
                    this.setF8_C_Cl(descF[i]);
            }
        }

        for (int d=0; d<AtomCouplesCC.length; d++) {

            int[] descF = new int[10];
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouplesNCl[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouplesNCl[d][1])) {

                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= 10) {
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouplesNCl[d][0].compareTo(AtomCouplesNCl[d][1]) == 0) {
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }

            for (int i=0; i<descF.length; i++) {
                if(i == 2)
                    this.setF3_N_Cl(descF[i]);
            }
        }




    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        this.setMLogP(0);
        try {
            DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
            this.setMLogP(descriptorMLogP.getMLogP());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setMLogP(MISSING_VALUE);
        }
    }

    private void CalculateACF(InsilicoMolecule Mol){
        this.setC_008(0);
        this.setO_060(0);
        try {
            DescriptorBlock acf = new ACF();
            acf.Calculate(Mol);
            this.setO_060(acf.GetByName("O-060").getValue());
            this.setC_008(acf.GetByName("C-008").getValue());
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
            this.setC_008(MISSING_VALUE);
            this.setO_060(MISSING_VALUE);
        }

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
}
