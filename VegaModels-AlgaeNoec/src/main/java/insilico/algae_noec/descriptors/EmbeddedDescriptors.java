package insilico.algae_noec.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.algae_noec.descriptors.weights.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.FunctionalGroups;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.Manipulator;
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

    public double Mw = MISSING_VALUE;
    public double MATS1p = MISSING_VALUE;
    public double MATS1s = MISSING_VALUE;
    public double GATS2v = MISSING_VALUE;
    public double MATS2v = MISSING_VALUE;
    public double GATS3s = MISSING_VALUE;
    public double ATSC3s = MISSING_VALUE;
    public double MATS4v = MISSING_VALUE;
    public double ALogP = MISSING_VALUE;
    public double MLogP = MISSING_VALUE;
    public double X2v = MISSING_VALUE;
    public double X1v = MISSING_VALUE;
    public double BEH3m = MISSING_VALUE;
    public double BEH6p = MISSING_VALUE;
    public double SpPosA_v = MISSING_VALUE;
    public double SpMAD_m = MISSING_VALUE;
    public double P_VSA_i_2 = MISSING_VALUE;
    public double piPC10 = MISSING_VALUE;
    public double PW5 = MISSING_VALUE;
    public double BIC1 = MISSING_VALUE;
    public double nCp = MISSING_VALUE;
    public double H047 = MISSING_VALUE;
    public double TCN = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
    }


    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateMw(Mol);

        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        MLogP = descriptorMLogP.MLogP;

        CalculateALogP(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateXv(Mol);
        CalculateBEH(Mol);
        CalculateP_VSA_i_2(Mol);
        CalculateWAP(Mol);
        CalculateBIC(Mol);
        CalculatenCp(Mol);
        CalculateACF(Mol);
        CalculateTCN(Mol);

    }

    private void CalculateTCN(InsilicoMolecule Mol){
        String[][] AtomCouples = {
                {"C", "N"}
        };

        TCN = 0;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            TCN = -999;
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            TCN = -999;
            return;
        }

        for (String[] atomCouple : AtomCouples) {
            for (int i = 0; i < nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(atomCouple[1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2)
                                TCN += TopoMat[i][j];
                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (atomCouple[0].compareTo(atomCouple[1]) == 0) {
                TCN /= 2;
            }
        }


    }

    private void CalculateACF(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorBlock descriptorBlock = new ACF();
        descriptorBlock.Calculate(Mol);
        H047 = descriptorBlock.GetByName("H-047").getValue();
    }

    private void CalculatenCp(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorBlock descriptorBlock = new FunctionalGroups();
        descriptorBlock.Calculate(Mol);
        nCp = descriptorBlock.GetByName("nCp").getValue();
    }

    private void CalculateBIC(InsilicoMolecule Mol) {
        BIC1 = 0;
        int MaxLag = 1;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BIC1 = MISSING_VALUE;
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
            BIC1 = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);

        // Neighborhood indices

        double bic_denom = 0;
        for (IBond bnd : m.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag=1; CurLag<=MaxLag; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK);
            for (int i=0; i<nSK; i++) {
                IAtom atStart = m.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (TopoDistMat[i][j] == CurLag) {
                        IAtom atEnd =  m.getAtom(j);
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


            BIC1 = bic_denom==0 ? 0 : ic / bic_denom;
        }
    }

    private void CalculateWAP(InsilicoMolecule Mol) {
        int pi_path = 10;
        int pw_path = 5;

        piPC10 = 0;
        PW5 = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            piPC10 = MISSING_VALUE;
            PW5 = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            piPC10 = MISSING_VALUE;
            PW5 = MISSING_VALUE;
            return;
        }

        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();

        double  CurMPC=0;
        int[] AtomWalk = GetAtomsWalks(pi_path, AdjMatDbl);
        double[][] AtomPath = GetAtomsPaths(pi_path, m);

        for (int i=0; i<nSK; i++) {
            CurMPC += AtomPath[i][1];
        }

        CurMPC /= 2.0;
        CurMPC = Math.log(1+CurMPC);

        piPC10 = CurMPC;

        double CurPW=0;
        AtomWalk = GetAtomsWalks(pw_path, AdjMatDbl);
        AtomPath = GetAtomsPaths(pw_path, m);

        for (int i=0; i<nSK; i++) {

            CurPW += (double)AtomPath[i][0] / (double)AtomWalk[i];
        }



        CurPW /= nSK;

        PW5 = CurPW;



    }

    private void CalculateP_VSA_i_2(InsilicoMolecule Mol) {

        P_VSA_i_2 = 0;
        short WEIGHT_I_IDX = 2; short bin_size = 2;

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(e.getMessage());
            P_VSA_i_2 = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();


        // Calculate VSA
        double[] VSA = new double[nSK];

        for (int i=0; i<nSK; i++) {

            IAtom at = m.getAtom(i);

            // R - van der waals radius
            double vdwR = GetVdWRadius(m, at);
            if (vdwR == MISSING_VALUE) continue;

            double coef = 0;
            for (IAtom connAt : m.getConnectedAtomsList(at)) {

                double connAt_vdwR = GetVdWRadius(m, connAt);
                if (connAt_vdwR == MISSING_VALUE) continue;

                // refR - reference bond length
                double refR = this.GetRefBondLength(at, connAt);
                if (refR == MISSING_VALUE) continue;

                // c - correction
                double c = 0;
                double bnd = MoleculeUtilities.Bond2Double(m.getBond(at, connAt));
                if (bnd == 1.5) c = 0.1;
                if (bnd == 2) c = 0.2;
                if (bnd == 3) c = 0.3;

                double g_1 = Math.max(  Math.abs(vdwR - connAt_vdwR) , GetRefBondLength(at, connAt) - c ) ;
                double g_2 = vdwR + connAt_vdwR;
                double g  = Math.min(g_1, g_2) ;

                coef += ( Math.pow(connAt_vdwR,2) - Math.pow( (vdwR - g), 2) ) / (g) ;
            }

            VSA[i] = ( 4.0 * Math.PI * Math.pow(vdwR,2) ) - Math.PI * vdwR * coef ;
        }



        // Cycle for all found weighting schemes

        // Sets needed weights
        double[] w = IonizationPotential.getWeightsNormalized(m);

        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE) continue;
            if ( bin_size == CalculateBin(w[i]))
                P_VSA_i_2 += VSA[i];
        }
    }

    private void CalculateBEH(InsilicoMolecule Mol) {
        BEH3m = 0; BEH6p = 0;

        int MAXEIG = 8;

        IAtomContainer m;

        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            BEH3m = MISSING_VALUE;
            BEH6p = MISSING_VALUE;
            return;
        }

        double[][] BurdenMat;
        BurdenMat = BurdenMatrix.getMatrix(m);


        int nSK = m.getAtomCount();
        double[] w;
        w = Mass.getWeights(m);

        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        double[] BEH = new double[MAXEIG];
        for (int i=0; i<MAXEIG; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
                valL = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
                valL = eigenvalues[i];
            }
            BEH[i] = valH;
        }

        BEH3m = BEH[2];

        // SpMAD
        double meanEig = 0;
        for (double e : eigenvalues) {
            meanEig += e;
        }
        meanEig /= eigenvalues.length;

        double SpAD = 0;
        for (double e : eigenvalues)
            SpAD += Math.abs(e - meanEig);

        SpMAD_m = SpAD / eigenvalues.length;


        w = Polarizability.getWeights(m);

        MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        BEH = new double[MAXEIG];
        for (int i=0; i<MAXEIG; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
                valL = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
                valL = eigenvalues[i];
            }
            BEH[i] = valH;
        }

        BEH6p = BEH[5];

        w = VanDerWaals.getWeights(m);

        MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;


        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);



        // SpAD, SpMAD, SpPosA
        meanEig = 0;
        double posSum = 0;
        for (double e : eigenvalues) {
            meanEig += e;
            if (e>0) posSum += e;
        }
        meanEig /= eigenvalues.length;

        SpAD = 0;
        for (double e : eigenvalues)
            SpAD += Math.abs(e - meanEig);


        SpPosA_v = posSum / nSK;
    }

    private void CalculateXv(InsilicoMolecule Mol) {
        X1v = 0; X2v = 0;
        int MaxPath1 = 1;
        int MaxPath2 = 2;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            X1v = -999;
            X2v = -999;
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            X1v = -999;
            X2v = -999;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);
        double[] curDescX = new double[MaxPath1+1];
        double[] curDescXv = new double[MaxPath1+1];
        double[] curDescXsol = new double[MaxPath1+1];

        // checks for missing weights
        for (int qnumber : Qnumbers)
            if (qnumber == -999)
                return;
        for (double v : ValenceVD)
            if (v == -999)
                return;

        // clears VD matrix from linked F
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ((ConnAugMatrix[i][j]>0) && (ConnAugMatrix[j][j]==9))
                    VD[i]--;
            }


        for (int k=0; k<MaxPath1; k++) {
            curDescX[k] = 0;
            curDescXv[k] = 0;
            curDescXsol[k] = 0;
        }

        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account

            // path 0
            curDescX[0] += Math.pow(VD[i], -0.5);
            curDescXv[0] += Math.pow(ValenceVD[i], -0.5);
            curDescXsol[0] += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);

            // path 1 - MaxPath
            for (int path=1; path<(MaxPath1+1); path++) {

                if (curDescX[path] == -999) continue;

                IAtom at = m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (List<IAtom> curPath : CurPaths) {
                    double prodXv = 1;
                    for (IAtom iAtom : curPath) {
                        int atIdx = m.indexOf(iAtom);
                        prodXv *= ValenceVD[atIdx];
                    }
                    curDescXv[path] += Math.pow(prodXv, -0.5);

                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath1+1); i++) {
            curDescXv[i] /= 2;
        }

        X1v = curDescXv[1];


        // 2

        curDescX = new double[MaxPath2+1];
        curDescXv = new double[MaxPath2+1];
        curDescXsol = new double[MaxPath2+1];

        // checks for missing weights
        for (int qnumber : Qnumbers)
            if (qnumber == -999)
                return;
        for (double v : ValenceVD)
            if (v == -999)
                return;

        // clears VD matrix from linked F
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ((ConnAugMatrix[i][j]>0) && (ConnAugMatrix[j][j]==9))
                    VD[i]--;
            }


        for (int k=0; k<MaxPath2; k++) {
            curDescX[k] = 0;
            curDescXv[k] = 0;
            curDescXsol[k] = 0;
        }

        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account

            // path 0
            curDescX[0] += Math.pow(VD[i], -0.5);
            curDescXv[0] += Math.pow(ValenceVD[i], -0.5);
            curDescXsol[0] += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);

            // path 1 - MaxPath
            for (int path=1; path<(MaxPath2+1); path++) {

                if (curDescX[path] == -999) continue;

                IAtom at = m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (List<IAtom> curPath : CurPaths) {
                    double prodXv = 1;
                    for (IAtom iAtom : curPath) {
                        int atIdx = m.indexOf(iAtom);
                        prodXv *= ValenceVD[atIdx];
                    }
                    curDescXv[path] += Math.pow(prodXv, -0.5);

                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath2+1); i++) {
            curDescXv[i] /= 2;
        }

        X2v = curDescXv[2];

    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol) {
        MATS1p = 0;
        MATS1s = 0;
        GATS2v = 0;
        MATS2v = 0;
        GATS3s = 0;
        ATSC3s = 0;
        MATS4v = 0;

        int PARAMETER_LAG_01 = 1;
        int PARAMETER_LAG_02 = 2;
        int PARAMETER_LAG_03 = 3;
        int PARAMETER_LAG_04 = 4;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MATS1p = MISSING_VALUE;
            MATS1s = MISSING_VALUE;
            GATS2v = MISSING_VALUE;
            MATS2v = MISSING_VALUE;
            GATS3s = MISSING_VALUE;
            ATSC3s = MISSING_VALUE;
            MATS4v = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS1p = MISSING_VALUE;
            MATS1s = MISSING_VALUE;
            GATS2v = MISSING_VALUE;
            MATS2v = MISSING_VALUE;
            GATS3s = MISSING_VALUE;
            ATSC3s = MISSING_VALUE;
            MATS4v = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w;

        // s
        try {
            EState ES = new EState(Mol.GetStructure());
            w = ES.getIS();
        } catch (Exception e) {
            w = new double[nSK];
            for (int i = 0; i < nSK; i++) w[i] = Descriptor.MISSING_VALUE;
        }
        boolean MissingWeight = false;
        for (int i = 0; i < nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;


        // Calculates weights averages
        double wA = 0;
        for (int i = 0; i < nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double ACS = 0, GearyAC = 0;
        double denom = 0, delta = 0;

        for (int i = 0; i < nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j = 0; j < nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_03) {
                    ACS += Math.abs((w[i] - wA) * (w[j] - wA));
                    GearyAC += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                GearyAC = 0;
            } else {
                GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double) (nSK - 1))) * denom);
            }
        }

        ACS /= 2.0;
        GATS3s = GearyAC;
        ATSC3s = ACS;

        // Calculates autocorrelations

        double MoranAC = 0;
        denom = 0; delta = 0;

        for (int i = 0; i < nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j = 0; j < nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_01) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double) nSK)) * denom);
            }
        }



        MATS1s = MoranAC;

        // V
        w = VanDerWaals.getWeights(m);
        MissingWeight = false;
        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
        }

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        MoranAC=0; GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_02) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    GearyAC += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MoranAC = 1;
                GearyAC = 0;
            } else {
                MoranAC = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
                GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
            }
        }

        GATS2v = GearyAC;
        MATS2v = MoranAC;

        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_04) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    GearyAC += Math.pow((w[i] - w[j]), 2);
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

        MATS4v = MoranAC;

        // p
        w = Polarizability.getWeights(m);
        MissingWeight = false;
        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
        }

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG_01) {
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

        MATS1p = MoranAC;
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        ALogP = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ALogP = MISSING_VALUE;
            log.warn(e.getMessage());
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);

        double[] Frags = acf.GetAllValues();

        for (double d : Frags)
            if (d == MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            ALogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }
    }

    private void CalculateMw(InsilicoMolecule Mol) {

        Mw = 0;
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            Mw = -999;
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[] H = new int[nSK];

            for (int i = 0; i < nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) {
                    log.warn(e.getMessage());
                }
            }

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");

            for (int i = 0; i < nSK; i++) {
                if (wMass[i] == -999)
                    Mw = -999;
            }

            for (int i = 0; i < nSK; i++) {
                if (Mw != -999) {
                    Mw += wMass[i];
                    if (H[i] > 0) {
                        Mw += HMass * H[i];
                    }
                }
            }

        } catch (Throwable e) {
            log.warn(e.getMessage());
            Mw = -999;
        }

    }

    private double GetVdWRadius(IAtomContainer m, IAtom at) {

        String s = at.getSymbol();
        if (s.equalsIgnoreCase("C"))
            return 1.950;
        if (s.equalsIgnoreCase("N"))
            return 1.950;
        if (s.equalsIgnoreCase("F"))
            return  1.496;
        if (s.equalsIgnoreCase("P"))
            return  2.287;
        if (s.equalsIgnoreCase("S"))
            return  2.185;
        if (s.equalsIgnoreCase("Cl"))
            return  2.044;
        if (s.equalsIgnoreCase("Br"))
            return  2.166;
        if (s.equalsIgnoreCase("I"))
            return 2.358;

        if (s.equalsIgnoreCase("O")) {

            // oxide
            if (m.getConnectedAtomsList(at).size() == 1)
                return 1.810;

            // acid
            if (Manipulator.CountImplicitHydrogens(at) > 0)
                return 2.152;

            return 1.779;
        }

        if (s.equalsIgnoreCase("H")) {
            IAtom connAt = m.getConnectedAtomsList(at).get(0);
            if (connAt.getSymbol().equalsIgnoreCase("0"))
                return 0.8;
            if (connAt.getSymbol().equalsIgnoreCase("N"))
                return 0.7;
            if (connAt.getSymbol().equalsIgnoreCase("P"))
                return 0.7;
            return 1.485;
        }

        return Descriptor.MISSING_VALUE;
    }

    private double GetRefBondLength(IAtom at1, IAtom at2) {
        double len = Descriptor.MISSING_VALUE;
        for (int i=0; i<RefBondLengths.length; i++) {
            if (AtomCouple(at1, at2, (String)RefBondLengths[i][0], (String)RefBondLengths[i][1])) {
                len = (Double)RefBondLengths[i][2];
                break;
            }
        }
        return len;
    }

    private final static Object[][] RefBondLengths = {
            {"Br","Br",2.54},
            {"Br","C",1.97},
            {"Br","Cl",2.36},
            {"Br","F",1.85},
            {"Br","H",1.44},
            {"Br","I",2.65},
            {"Br","N",1.84},
            {"Br","O",1.58},
            {"Br","P",2.37},
            {"Br","S",2.21},
            {"C","C",1.54},
            {"C","Cl",1.8},
            {"C","F",1.35},
            {"C","H",1.06},
            {"C","I",2.12},
            {"C","N",1.47},
            {"C","O",1.43},
            {"C","P",1.85},
            {"C","S",1.81},
            {"Cl","Cl",2.31},
            {"Cl","F",1.63},
            {"Cl","H",1.22},
            {"Cl","I",2.56},
            {"Cl","N",1.74},
            {"Cl","O",1.41},
            {"Cl","P",2.01},
            {"Cl","S",2.07},
            {"F","F",1.28},
            {"F","H",0.87},
            {"F","I",2.04},
            {"F","N",1.41},
            {"F","O",1.32},
            {"F","P",1.5},
            {"F","S",1.64},
            {"H","I",1.63},
            {"H","N",1.01},
            {"H","O",0.97},
            {"H","P",1.41},
            {"H","S",1.31},
            {"I","I",2.92},
            {"I","N",2.26},
            {"I","O",2.14},
            {"I","P",2.49},
            {"I","S",2.69},
            {"N","N",1.45},
            {"N","O",1.46},
            {"N","P",1.6},
            {"N","S",1.76},
            {"O","O",1.47},
            {"O","P",1.57},
            {"O","S",1.57},
            {"P","P",2.26},
            {"P","S",2.07},
            {"S","S",2.05}
    };

    private int CalculateBin(double value) {

        // Ionization potential
        if (value < 1) return 1;
        if (value < 1.15) return 2;
        if (value < 1.25) return 3;
        return 4;


    }

    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
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
}
