package insilico.fish_lc50.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.fish_lc50.descriptors.weights.*;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Log4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;
    private double ATSC1s = MISSING_VALUE;
    private double MATS5s = MISSING_VALUE;
    private double BEH3m = MISSING_VALUE;
    private double BEH4m = MISSING_VALUE;
    private double BEH7p = MISSING_VALUE;
    private double CATS2D_2_DA = MISSING_VALUE;
    private double X2sol = MISSING_VALUE;
    private double Mp = MISSING_VALUE;
    private double NPerc = MISSING_VALUE;
    private double BIC3 = MISSING_VALUE;
    private double ALogP = MISSING_VALUE;
    private double MLogP = MISSING_VALUE;

    private double Mw = MISSING_VALUE;

    private final static String TYPE_A = "A";
    private final static String TYPE_D = "D";
    private final static String[][] AtomCouples = {
            {TYPE_D, TYPE_A},
    };
    private double[][] ConnAugMatrixCATS;
    

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateAutocorrelation(Mol);
        CalculateBE(Mol);
        CalculateCATS2D(Mol);
        CalculateCI(Mol);
        CalculateMLogP(Mol);
        CalculateALogP(Mol);
        CalculateConstitutional(Mol);
        CalculateIC(Mol);
    }

    private void CalculateConstitutional(InsilicoMolecule Mol){
        this.setMw(0);
        this.setNPerc(0);
        this.setMp(0);

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMw(MISSING_VALUE);
            this.setNPerc(MISSING_VALUE);
            this.setMp(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nTotH=0;
            int nC=0, nN=0, nO=0, nP=0, nS=0;
            int nI=0, nF=0, nCl=0, nBr=0, nB=0;
            int nHet=0;
            double mw=0, amw=0, sv=0, mv=0, sp=0, mp=0, se=0, me=0;


            //// Counts on atoms

            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];


                if (CurAt.getSymbol().equalsIgnoreCase("C"))
                    nC++;
                else
                    nHet++;

                if (CurAt.getSymbol().equalsIgnoreCase("N"))
                    nN++;
                if (CurAt.getSymbol().equalsIgnoreCase("O"))
                    nO++;
                if (CurAt.getSymbol().equalsIgnoreCase("P"))
                    nP++;
                if (CurAt.getSymbol().equalsIgnoreCase("S"))
                    nS++;
                if (CurAt.getSymbol().equalsIgnoreCase("F"))
                    nF++;
                if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                    nCl++;
                if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                    nBr++;
                if (CurAt.getSymbol().equalsIgnoreCase("I"))
                    nI++;
                if (CurAt.getSymbol().equalsIgnoreCase("B"))
                    nB++;

            }

            this.setNPerc((nN/(double)(nSK + nTotH))*100);

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");
            double[] wVdW = VanDerWaals.getWeights(curMol);
            double HVdW = VanDerWaals.GetVdWVolume("H");
            double[] wPol = Polarizability.getWeights(curMol);
            double HPol = Polarizability.GetPolarizability("H");
            double[] wEl = Electronegativity.getWeights(curMol);
            double HEl = Electronegativity.GetElectronegativity("H");



            for (int i=0; i<nSK; i++) {
                if (wMass[i] == -999)
                    mw = -999;
                if (wVdW[i] == -999)
                    sv = -999;
                if (wPol[i] == -999)
                    sp = -999;
                if (wEl[i] == -999)
                    se = -999;
            }

            for (int i=0; i<nSK; i++) {
                if (mw != -999) {
                    mw += wMass[i];
                    if (H[i]>0) {
                        mw += HMass * H[i];
                    }
                }
                if (sv != -999) {
                    sv += wVdW[i];
                    if (H[i]>0)
                        sv += HVdW * H[i];
                }
                if (sp != -999) {
                    sp += wPol[i];
                    if (H[i]>0)
                        sp += HPol * H[i];
                }
                if (se != -999) {
                    se += wEl[i];
                    if (H[i]>0)
                        se += HEl * H[i];
                }
            }

            if (sp != -999)
                mp = sp/(nSK + nTotH);

            this.setMw(mw);
            this.setMp(mp);
        } catch (Throwable e) {
            this.setMw(MISSING_VALUE);
            this.setNPerc(MISSING_VALUE);
            this.setMp(MISSING_VALUE);
        }
    }

    private void CalculateIC(InsilicoMolecule Mol){
        this.setBIC3(0);

        int MaxPath = 4;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setBIC3(MISSING_VALUE);
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
            setBIC3(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);


        // Information content

        int[] TopoDistFreq = new int[nSK];  // frequencies of topological distances
        double TopoDistFreqSum = 0;
        for (int i=0; i<nSK; i++)
            TopoDistFreq[i] = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i+1; j<nSK; j++) {
                TopoDistFreq[TopoDistMat[i][j]]++;
                TopoDistFreqSum += TopoDistMat[i][j];
            }

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

                        List<IAtom> sp = getShortestPath(m, atStart, atEnd);

                        StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                        for (int k=0; k<(sp.size()-1); k++) {
                            int a = m.indexOf(sp.get(k));
                            int b = m.indexOf(sp.get(k + 1));
                            if (ConnMat[a][b] == 1)
                                bufPath.append("s");
                            if (ConnMat[a][b] == 2)
                                bufPath.append("d");
                            if (ConnMat[a][b] == 3)
                                bufPath.append("t");
                            if (ConnMat[a][b] == 1.5)
                                bufPath.append("a");
                            bufPath.append(sp.get(k + 1).getAtomicNumber());
                            bufPath.append("(").append(VertexDeg[m.indexOf(sp.get(k + 1))]).append(")");
                        }
                        CurNeig.add(bufPath.toString());
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

            double cic = (Math.log(nSK) / Math.log(2)) - ic;

            double bic = bic_denom==0 ? 0 : ic / bic_denom;


            if(CurLag == 3)
                this.setBIC3(bic);
        }


    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        this.setMLogP(descriptorMLogP.MLogP);

    }

    private void CalculateALogP(InsilicoMolecule Mol){
        this.setALogP(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setALogP(MISSING_VALUE);
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF = new ACF();
        CurACF.Calculate(Mol);

        double LogP = 0;
        double[] Frags = CurACF.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }

        this.setALogP(LogP);
    }



    private void CalculateCI(InsilicoMolecule Mol){

        this.setX2sol(0);

        int MaxPath = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setX2sol(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setX2sol(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);
        double[] curDescX = new double[MaxPath+1];
        double[] curDescXv = new double[MaxPath+1];
        double[] curDescXsol = new double[MaxPath+1];

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


        for (int k=0; k<MaxPath; k++) {
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
            for (int path=1; path<(MaxPath+1); path++) {

                if (curDescX[path] == -999) continue;

                IAtom at = m.getAtom(i);
                List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(m, at, path);
                for (List<IAtom> curPath : CurPaths) {
                    double prodX = 1;
                    double prodXv = 1;
                    int prodQuantum = 1;
                    for (IAtom iAtom : curPath) {
                        int atIdx = m.indexOf(iAtom);
//                        if (ConnMatrix[atIdx][atIdx] == 9)
//                            continue; // F not taken into account
                        prodX *= VD[atIdx];
                        prodXv *= ValenceVD[atIdx];
                        prodQuantum *= Qnumbers[atIdx];
                    }
                    curDescX[path] += Math.pow(prodX, -0.5);
                    curDescXv[path] += Math.pow(prodXv, -0.5);
                    curDescXsol[path] += (1.00 / Math.pow(2.00, (double) (path + 1))) *
                            ((double) prodQuantum) * Math.pow(prodX, -0.5);
                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath+1); i++) {
            curDescX[i] /= 2;
            curDescXv[i] /= 2;
            curDescXsol[i] /= 2;
        }

        // Sets descriptors
        for (int i=0; i<(MaxPath+1); i++) {
            if(i == 2)
                this.setX2sol(curDescXsol[i]);
        }
    }
    
    private void CalculateCATS2D(InsilicoMolecule Mol){
        this.setCATS2D_2_DA(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setCATS2D_2_DA(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_2_DA(MISSING_VALUE);
            return;
        }
        try {
            ConnAugMatrixCATS = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_2_DA(MISSING_VALUE);
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (String[] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[3];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], atomCouple[0])) {
                    for (int j = 0; j < nSK; j++) {
                        if (i == j) continue;
                        if (isIn(CatsTypes[j], atomCouple[1])) {

                            if (TopoMat[i][j] < 3)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++)
                if (i == 2)
                    this.setCATS2D_2_DA(desc[i]);

        }

    }

    private void CalculateBE(InsilicoMolecule Mol){
        setBEH3m(0); setBEH4m(0); setBEH7p(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setBEH3m(MISSING_VALUE); setBEH4m(MISSING_VALUE); setBEH7p(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            setBEH3m(MISSING_VALUE); setBEH4m(MISSING_VALUE); setBEH7p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        // m
        double[] w = Mass.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
            }
            if(i==2)
                this.setBEH3m(valH);
            if(i==3)
                this.setBEH4m(valH);
        }

        w = Polarizability.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
            }
            if(i==6)
                this.setBEH7p(valH);
        }
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol){
        setATSC1s(0); setMATS5s(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setATSC1s(MISSING_VALUE); setMATS5s(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setATSC1s(MISSING_VALUE); setMATS5s(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        double[] w;

        try {
            EState ES = new EState(Mol.GetStructure());
            w = ES.getIS();
        } catch (Exception e) {
            w = new double[nSK];
            for (int i=0; i<nSK; i++) w[i]=Descriptor.MISSING_VALUE;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double AC=0, ACS=0;

        for (int i=0; i<nSK; i++) {

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 1) {
                    ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                }
        }


        ACS /= 2.0;

        // Sets descriptors
        this.setATSC1s(ACS);

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 5) {
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

        setMATS5s(MoranAC);
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
                if (ConnAugMatrixCATS[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrixCATS[i][j]==1) {
                            if (ConnAugMatrixCATS[j][j] == 8) {
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1)
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrixCATS[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrixCATS[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrixCATS[i][j] == 2) && (ConnAugMatrixCATS[j][j] == 8) )
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] > 1.5) ) {
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] != 1) ) {
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
            if (tA) AtomTypes[i].add(TYPE_A);
            if (tD) AtomTypes[i].add(TYPE_D);

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



}
