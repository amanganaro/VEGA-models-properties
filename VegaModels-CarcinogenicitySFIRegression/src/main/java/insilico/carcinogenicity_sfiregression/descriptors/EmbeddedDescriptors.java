package insilico.carcinogenicity_sfiregression.descriptors;

import Jama.Matrix;
import insilico.carcinogenicity_sfiregression.descriptors.weights.ACF;
import insilico.carcinogenicity_sfiregression.descriptors.weights.FunctionalGroups;
import insilico.carcinogenicity_sfiregression.descriptors.weights.Mass;
import insilico.carcinogenicity_sfiregression.descriptors.weights.Polarizability;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.weight.VertexDegree;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Slf4j
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;
    public double C_041 = MISSING_VALUE;
    public double ATS8m = MISSING_VALUE;
    public double GATS6p = MISSING_VALUE;
    public double CATS2D_3_DL = MISSING_VALUE;
    public double CATS2D_7_DL = MISSING_VALUE;
    public double nN_N = MISSING_VALUE;
    public double IC4 = MISSING_VALUE;
    public double B2Cl_Cl = MISSING_VALUE;
    public double B4O_Cl = MISSING_VALUE;
    public double B7Cl_Cl = MISSING_VALUE;
    public double B8Cl_Cl = MISSING_VALUE;
    public double SRW7 = MISSING_VALUE;


    private final static String TYPE_L = "L";
    private final static String TYPE_D = "D";

    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_L},
    };

    private double[][] ConnAugMatrix;

    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateCats2D(Mol);
        CalculateFG(Mol);
        CalculateIC(Mol);
        CalculateTopoDistances(Mol);
        CalculateWAP(Mol);
    }

    private void CalculateWAP(InsilicoMolecule Mol){
        SRW7 = 0;
        short path = 7;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SRW7 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SRW7 = MISSING_VALUE;
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];


//        if (getBoolProperty(PARAMETER_INCLUDE_SRW)) {
//            int SRWMaxPath = pathsList.get(pathsList.size()-1);
        int[] MolSRW = GetSRWToLag(path, AdjMatDbl);
        for (int i=1; i<=path; i++)
            if(i == 7)
                SRW7 = MolSRW[i];
//        }
    }

    private void CalculateTopoDistances(InsilicoMolecule Mol){

        B2Cl_Cl = 0; B4O_Cl = 0; B7Cl_Cl = 0; B8Cl_Cl = 0;

        String AtomCouples1[][] = {
                {"O", "Cl"},
        };

        String AtomCouples2[][] = {
                {"Cl", "Cl"},
        };

        int MAX_TOPO_DISTANCE = 10;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            B2Cl_Cl = MISSING_VALUE; B4O_Cl = MISSING_VALUE;
            B7Cl_Cl = MISSING_VALUE; B8Cl_Cl = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            B2Cl_Cl = MISSING_VALUE; B4O_Cl = MISSING_VALUE;
            B7Cl_Cl = MISSING_VALUE; B8Cl_Cl = MISSING_VALUE;
            return;
        }

        for (int d=0; d<AtomCouples1.length; d++) {

            int descT = 0;
            int[] descB = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouples[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouples[d][1])) {

                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j]-1] = 1;
                            }

                        }
                    }
                }
            }

            for (int i=0; i<descB.length; i++)
            {
                if(i == 4)
                    B4O_Cl = descB[i];
            }


        }

        for (int d=0; d<AtomCouples2.length; d++) {

            int[] descB = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouples[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouples[d][1])) {

                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j]-1] = 1;
                            }

                        }
                    }
                }
            }

            for (int i=0; i<descB.length; i++)
            {
                if(i == 2)
                    B2Cl_Cl = descB[i];
                if(i == 7)
                    B7Cl_Cl = descB[i];
                if(i == 8)
                    B8Cl_Cl = descB[i];
            }

        }
    }

    private void CalculateIC(InsilicoMolecule Mol) {

        IC4 = 0;
        int MaxPath = 4;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            IC4 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] ConnMat;
        int[][] TopoDistMat;
        try {
            ConnMat = Mol.GetMatrixConnectionAugmented();
            TopoDistMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            IC4 = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);


        // Neighborhood indices

        double bic_denom = 0;
        for (IBond bnd : m.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag=1; CurLag<=MaxPath; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK);
            for (int i=0; i<nSK; i++) {
                IAtom atStart = m.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (TopoDistMat[i][j] == CurLag) {
                        IAtom atEnd =  m.getAtom(j);

                        // This should be the correct way to refactor the old code
                        // but with new CDK results are different compared to CDK 1.4.9
                        // so we use the old getShortestPath method, its source code has been taken
                        // and put as static method in this same class

//                        ShortestPaths sps = new ShortestPaths(m, atStart);
//                        List<IAtom> sp = Arrays.asList(sps.atomsTo(atEnd));

                        // OLD DEPRECATED METHOD
                        List<IAtom> sp = getShortestPath(m, atStart, atEnd);

                        String bufPath = "" + sp.get(0).getAtomicNumber();
                        for (int k=0; k<(sp.size()-1); k++) {
                            int a = m.indexOf(sp.get(k));
                            int b = m.indexOf(sp.get(k + 1));
                            if (ConnMat[a][b] == 1)
                                bufPath += "s";
                            if (ConnMat[a][b] == 2)
                                bufPath += "d";
                            if (ConnMat[a][b] == 3)
                                bufPath += "t";
                            if (ConnMat[a][b] == 1.5)
                                bufPath += "a";
                            bufPath += sp.get(k+1).getAtomicNumber();
                            bufPath += "(" + VertexDeg[m.indexOf(sp.get(k + 1))] + ")";
                        }
                        CurNeig.add(bufPath);
                    }
                }
                Collections.sort(CurNeig);
                NeigList.add(CurNeig);
            }

            // Calculates equivalence classes
            ArrayList<ArrayList<String>> G = new ArrayList<>();
            ArrayList<Integer> Gn = new ArrayList<>();
            for (int i=0; i<nSK; i++) {
                ArrayList<String> CurNeig = NeigList.get(i);
                boolean foundMatch = false;
                for (int k=0; k<G.size(); k++) {
                    if (CompareNeigVector(CurNeig, G.get(k))) {
                        foundMatch = true;
                        int buf = Gn.get(k);
                        Gn.set(k, (buf+1));
                        break;
                    }
                }
                if (!foundMatch) {
                    G.add(CurNeig);
                    Gn.add(1);
                }
            }

            // Calculate IC and CIC indices

            double ic=0;
            for (int i=0; i<Gn.size(); i++)
                ic += ((double)Gn.get(i)/nSK) * (Math.log((double)Gn.get(i)/nSK));
            ic = (-1.00 / Math.log(2)) * ic;

            if(CurLag == 4)
                IC4 = ic;

        }

    }

    private void CalculateFG(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        nN_N = 0;
        DescriptorBlock fg = new FunctionalGroups();
        fg.Calculate(Mol);
        nN_N = fg.GetByName("nN-N").getValue();
    }

    private void CalculateCats2D(InsilicoMolecule Mol) {
        CATS2D_3_DL = 0; CATS2D_7_DL = 0;
        int MAX_CATS_DISTANCE = 10;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            CATS2D_3_DL = MISSING_VALUE; CATS2D_7_DL = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DL = MISSING_VALUE; CATS2D_7_DL = MISSING_VALUE;
            return;
        }
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            CATS2D_3_DL = MISSING_VALUE; CATS2D_7_DL = MISSING_VALUE;
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (String[] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[MAX_CATS_DISTANCE];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (isIn(CatsTypes[j], atomCouple[1])) {

                            if (TopoMat[i][j] < MAX_CATS_DISTANCE)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++) {
                if (i == 3)
                    CATS2D_3_DL = desc[i];
                if (i == 7)
                    CATS2D_7_DL = desc[i];
            }
        }
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        ACF acf = new ACF();
        acf.Calculate(Mol);
        C_041 = acf.GetByName("C-041").getValue();
    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){
        ATS8m = 0; GATS6p = 0;
        double[] w_p, w_m;

        short lag_m = 8;
        short lag_p = 6;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ATS8m = MISSING_VALUE; GATS6p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            ATS8m = MISSING_VALUE; GATS6p = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        w_p = Polarizability.getWeights(m);
        w_m = Mass.getWeights(m);;


        double wA_p = 0;
        double wA_m = 0;
        for (int i=0; i<nSK; i++) {
            wA_p += w_p[i];
            wA_m += w_m[i];
        }

        wA_p = wA_p / ((double) nSK);
        wA_m = wA_m / ((double) nSK);

        double AC=0;
        double denom = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_m[i] - wA_m), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_m) {
                    AC += w_m[i] * w_m[j];
                }
        }

        AC /= 2.0;
        AC = Math.log(1 + AC);

        ATS8m = AC;

        double GearyAC=0;
        denom = 0;
        double delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w_p[i] - wA_p), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag_p) {
                    GearyAC += Math.pow((w_p[i] - w_p[j]), 2);
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

        GATS6p = GearyAC;
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
                if (ConnAugMatrix[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrix[i][j]==1) {
                            if (ConnAugMatrix[j][j] == 8) {
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
                    if (ConnAugMatrix[i][j]>0) {
                        if (ConnAugMatrix[i][j] == 1)
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
                    if (ConnAugMatrix[i][j]>0) {
                        if (ConnAugMatrix[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrix[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrix[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
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
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] > 1.5) ) {
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
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] != 1) ) {
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

    private boolean CompareNeigVector(ArrayList<String> A, ArrayList<String> B) {

        // note: integer in the vectors should be already sorted

        if (A.size()!=B.size())
            return false;
        else {
            for (int i=0; i<A.size(); i++)
                if (!A.get(i).equalsIgnoreCase(B.get(i)))
                    return false;
        }

        return true;
    }

    /**
     * This method is taken from CDK 1.4.9 and used for these descriptors for retro-compatibility
     *
     * @param atomContainer
     * @param start
     * @param end
     * @return
     */
    private static List<IAtom> getShortestPath(IAtomContainer atomContainer, IAtom start, IAtom end) {
        int natom = atomContainer.getAtomCount();
        int endNumber = atomContainer.indexOf(end);
//        int endNumber = atomContainer.getAtomNumber(end);
        int startNumber = atomContainer.indexOf(start);
//        int startNumber = atomContainer.getAtomNumber(start);
        int[] dist = new int[natom];
        int[] previous = new int[natom];
        for (int i = 0; i < natom; i++) {
            dist[i] = 99999999;
            previous[i] = -1;
        }
        dist[atomContainer.indexOf(start)] = 0;
//        dist[atomContainer.getAtomNumber(start)] = 0;

        List<IAtom> Slist = new ArrayList<IAtom>();
        List<Integer> Qlist = new ArrayList<Integer>();
        for (int i = 0; i < natom; i++) Qlist.add(i);

        while (true) {
            if (Qlist.size() == 0) break;

            // extract min
            int u = 999999;
            int index = 0;
            for (Integer tmp : Qlist) {
                if (dist[tmp] < u) {
                    u = dist[tmp];
                    index = tmp;
                }
            }
            Qlist.remove(Qlist.indexOf(index));
            Slist.add(atomContainer.getAtom(index));
            if (index == endNumber) break;

            // relaxation
            List<IAtom> connected = atomContainer.getConnectedAtomsList(atomContainer.getAtom(index));
            for (IAtom aConnected : connected) {
                int anum = atomContainer.indexOf(aConnected);
//                int anum = atomContainer.getAtomNumber(aConnected);
                if (dist[anum] > dist[index] + 1) { // all edges have equals weights
                    dist[anum] = dist[index] + 1;
                    previous[anum] = index;
                }
            }
        }

        ArrayList<IAtom> tmp = new ArrayList<IAtom>();
        int tmpSerial = endNumber;
        while (true) {
            tmp.add(0, atomContainer.getAtom(tmpSerial));
            tmpSerial = previous[tmpSerial];
            if (tmpSerial == startNumber) {
                tmp.add(0, atomContainer.getAtom(tmpSerial));
                break;
            }
        }
        return tmp;
    }



    private static double Log(int base,double x) {
        double Logbx = Math.log10(x)/Math.log10((double)base);
        return Logbx;
    }

    private int[] GetAtomsWalks(int WalksOrder, double[][] AdjMatrix) {

        int[] walks = new int[AdjMatrix.length];

        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);

        for (int k=1; k<WalksOrder; k++) {
            mWalks = mWalks.times(mAdj);
        }

        for (int i=0; i<AdjMatrix.length; i++) {
            int CurSum = 0;
            for (int j=0; j<AdjMatrix[0].length; j++)
                CurSum += mWalks.get(i, j);
            walks[i] = CurSum;
        }

        return walks;
    }

    // index 0 = count of atom paths
    // index 1 = sum of paths weighted by bond order
    private double[][] GetAtomsPaths(int PathsOrder, IAtomContainer Mol) {

        int nSK = Mol.getAtomCount();
        double[][] paths = new double[nSK][2];

        for (int i=0; i<nSK; i++) {
            IAtom a = Mol.getAtom(i);
            int PathNum = 0;
            double totBO = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j)
                    continue;
                IAtom b =  Mol.getAtom(j);
                List<List<IAtom>> list= PathTools.getAllPaths(Mol, a, b);
                for (int k=0; k<list.size(); k++) {
                    if (list.get(k).size() == (PathsOrder+1)) {
                        PathNum++;

                        double curBO=1;
                        for (int at_idx=0; at_idx<list.get(k).size()-1; at_idx++) {
                            IAtom at_1 = list.get(k).get(at_idx);
                            IAtom at_2 = list.get(k).get(at_idx+1);
                            IBond curBond = Mol.getBond(at_2, at_1);
                            if (curBond.getFlag(CDKConstants.ISAROMATIC))
                                curBO *= 1.5;
                            else {
                                if (curBond.getOrder() == IBond.Order.SINGLE) curBO *= 1;
                                if (curBond.getOrder() == IBond.Order.DOUBLE) curBO *= 2;
                                if (curBond.getOrder() == IBond.Order.TRIPLE) curBO *= 3;
                                if (curBond.getOrder() == IBond.Order.QUADRUPLE) curBO *= 4;
                            }
                        }
                        totBO+=curBO;
                    }
                }
                paths[i][0] = PathNum;
                paths[i][1] = totBO;
            }
        }

        return paths;
    }


    private int[] GetSRWToLag(int Lag, double[][] AdjMatrix) {

        int nSK = AdjMatrix.length;
        int[] MolSRW = new int[Lag+1];

        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);

        for (int k=1; k<(Lag+1); k++) {

            MolSRW[k] = 0;
            for (int i=0; i<nSK; i++) {
                MolSRW[k] += mWalks.get(i, i);
            }

            mWalks = mWalks.times(mAdj);
        }

        return MolSRW;
    }



}
