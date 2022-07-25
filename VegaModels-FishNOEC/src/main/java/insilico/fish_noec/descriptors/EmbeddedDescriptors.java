package insilico.fish_noec.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.BurdenMatrix;
import insilico.fish_noec.descriptors.weights.*;
import insilico.fish_noec.descriptors.weights.ACF;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;

@Log4j
@Data
public class EmbeddedDescriptors {

    public final int MISSING_VALUE = -999;

    private double ATS5p = MISSING_VALUE;
    private double ALogP = MISSING_VALUE;
    private double X3sol = MISSING_VALUE;
    private double X1v = MISSING_VALUE;
    private double BEH4m = MISSING_VALUE;
    private double MW = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol) {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) {
        CalculateAutoCorrelation(Mol);
        CalculateALogP(Mol);
        CalculateConnectivityIndices(Mol);
        CalculateBurdenEigenValue(Mol);
        CalculateMW(Mol);
    }

    private void CalculateMW(InsilicoMolecule Mol){
        this.setMW(0);
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMW(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];


            double mw=0, amw=0, sv=0, mv=0, sp=0, mp=0, se=0, me=0;

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
            }

            this.setMW(mw);


        } catch (Throwable e) {
            this.setMW(MISSING_VALUE);
        }
    }

    private void CalculateBurdenEigenValue(InsilicoMolecule Mol){
         setBEH4m(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setBEH4m(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            setBEH4m(MISSING_VALUE);
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
            if(i==3)
                this.setBEH4m(valH);
        }
    }

    private void CalculateConnectivityIndices(InsilicoMolecule Mol){
        this.setX3sol(0); this.setX1v(0);

        int MaxPath = 5;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setX3sol(MISSING_VALUE); this.setX1v(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setX3sol(MISSING_VALUE); this.setX1v(MISSING_VALUE);
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
                    curDescXv[path] += Math.pow(prodXv, -0.5);
                    curDescXsol[path] += (1.00 / Math.pow(2.00, (double) (path + 1))) *
                            ((double) prodQuantum) * Math.pow(prodX, -0.5);
                }

            }
        }

        // descriptors with path>0 counted all paths twice
        for (int i=1; i<(MaxPath+1); i++) {
            curDescXv[i] /= 2;
            curDescXsol[i] /= 2;
        }

        // Sets descriptors
        for (int i=0; i<(MaxPath+1); i++) {
            if (i == 1)
                this.setX1v(curDescXv[i]);
            if (i == 3)
                this.setX3sol(curDescXsol[i]);
        }
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        setALogP(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setALogP(MISSING_VALUE);
            return;
        }

        DescriptorBlock acf = new ACF();
        acf.Calculate(Mol);

        double LogP = 0;
        double[] Frags = acf.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }

        this.setALogP(LogP);
    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){
        setATS5p(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            setATS5p(MISSING_VALUE);
            return;
        }

        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            setATS5p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        //p
        double[] w = Polarizability.getWeights(m);
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double AC=0;

        for (int i=0; i<nSK; i++) {
            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 5) {
                    AC += w[i] * w[j];
                }
        }

        // AC transformed in log form
        AC /= 2.0;
        AC = Math.log(1 + AC);

        this.setATS5p(AC);
    }










}
