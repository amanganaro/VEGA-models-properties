package insilico.rba_irfmn.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.rba_irfmn.descriptors.weights.*;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Arrays;
import java.util.List;

@Log4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double X2v = MISSING_VALUE;
    private double MATS6m = MISSING_VALUE;
    private double MATS8v = MISSING_VALUE;
    private double MATS5p = MISSING_VALUE;
    private double BEH2e = MISSING_VALUE;
    private double BEH1p = MISSING_VALUE;
    private double nArOH = MISSING_VALUE;
    private double MLogP = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateMLogP(Mol);
        CalculateX2v(Mol);
        CalculateFG(Mol);
        CalculateAutoCorrelation(Mol);
        CalculateBEH(Mol);
    }
    
    private void CalculateBEH(InsilicoMolecule Mol){
        this.setBEH1p(0); this.setBEH2e(0);

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            this.setBEH1p(MISSING_VALUE); this.setBEH1p(MISSING_VALUE);            
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            this.setBEH1p(MISSING_VALUE); this.setBEH1p(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        // p
        double[] w = Polarizability.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            if (i == 0) {
                double valH, valL;
                if (i>(eigenvalues.length-1)) {
                    this.setBEH1p(0);
                } else {
                    this.setBEH1p(eigenvalues[eigenvalues.length-1-i]);
                }
            }
        }


        // e
        w = Electronegativity.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }
        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            if (i == 1) {
                if (i>(eigenvalues.length-1)) {
                    this.setBEH2e(0);
                } else {
                    this.setBEH2e(eigenvalues[eigenvalues.length-1-i]);
                }
            }
        }
    }

    private void CalculateAutoCorrelation(InsilicoMolecule Mol){
        this.setMATS5p(0); this.setMATS6m(0); this.setMATS8v(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMATS5p(MISSING_VALUE); this.setMATS6m(MISSING_VALUE); this.setMATS8v(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setMATS5p(MISSING_VALUE); this.setMATS6m(MISSING_VALUE); this.setMATS8v(MISSING_VALUE);
            return;
        }

        // !!! in origine era usata la topological matrix del cdk

        int nSK = m.getAtomCount();
        
        // p
        double[] w = Polarizability.getWeights(m);

        double wA = 0;
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
        
        this.setMATS5p(MoranAC);
        
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

        this.setMATS8v(MoranAC);

        //m
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

        this.setMATS6m(MoranAC);
        
    }

    private void CalculateFG(InsilicoMolecule Mol){
        this.setNArOH(0);
        try {
            DescriptorBlock FG = new FG();
            FG.Calculate(Mol);
            this.setNArOH(FG.GetByName("nArOH").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setNArOH(MISSING_VALUE);
        }
    }

    private void CalculateX2v(InsilicoMolecule Mol){

        this.setX2v(0);
        int MaxPath2 = 2;


        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setX2v(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setX2v(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);

        double[] curDescX = new double[MaxPath2+1];
        double[] curDescXv = new double[MaxPath2+1];
        double[] curDescXsol = new double[MaxPath2+1];

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

        this.setX2v(curDescXv[2]);

    }





    private void CalculateMLogP(InsilicoMolecule Mol){
        try {
            this.setMLogP(0);
            DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
            this.setMLogP(descriptorMLogP.getMLogP());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setMLogP(MISSING_VALUE);
        }

    }

}
