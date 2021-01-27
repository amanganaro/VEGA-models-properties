package insilico.skin_caesar.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.skin_caesar.descriptors.weights.*;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double nN = MISSING_VALUE;
    private double GNar = MISSING_VALUE;
    private double X2v = MISSING_VALUE;
    private double EEig10ri = MISSING_VALUE;
    private double GGI8 = MISSING_VALUE;
    private double nCconj = MISSING_VALUE;
    private double O_058 = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateFG(Mol);
        CalculateACF(Mol);
        CalculateConstitutional(Mol);
        CalculateTopological(Mol);
        CalculateX2v(Mol);
        CalculateEA(Mol);
        CalculateGGI(Mol);
    }

    private void CalculateGGI(InsilicoMolecule Mol){
        this.setGGI8(0);
        int MaxPath = 8;

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setGGI8(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setGGI8(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[][] AdjMat = AdjacencyMatrix.getMatrix(curMol);
            int[] VertexDegrees = VertexDegree.getWeights(curMol, false);

            Matrix mAdj = new Matrix(AdjMat.length, AdjMat[0].length);
            for (int i=0; i<AdjMat.length; i++)
                for (int j=0; j<AdjMat[0].length; j++)
                    mAdj.set(i, j, (double)AdjMat[i][j]);

            Matrix mRecSqrDist = new Matrix(TopoMat.length, TopoMat[0].length);
            for (int i=0; i<TopoMat.length; i++)
                for (int j=0; j<TopoMat[0].length; j++) {
                    double CurVal = 0;
                    if (TopoMat[i][j] != 0)
                        CurVal = 1 / Math.pow((double)TopoMat[i][j], 2);
                    mRecSqrDist.set(i, j, CurVal);
                }

            Matrix multMat = mAdj.times(mRecSqrDist);

            double[][] CTMatrix = new double[nSK][nSK];
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++) {
                    if (i == j)
                        CTMatrix[i][j] = VertexDegrees[i];
                    else {
                        CTMatrix[i][j] = multMat.get(i, j) - multMat.get(j, i);
                    }
                }

            int[] PathCount = new int[MaxPath];
            double[] GGIval = new double[MaxPath];
            for (int i=0; i<MaxPath; i++) {
                PathCount[i] = 0;
                GGIval[i] = 0;
            }

            for (int i=0; i<nSK; i++)
                for (int j=i; j<nSK; j++) {

                    int CurPath = TopoMat[i][j];
                    if  ((CurPath>0) && (CurPath<=MaxPath)) {
                        PathCount[CurPath-1]++;
                        GGIval[CurPath-1] += Math.abs(CTMatrix[i][j]);
                    }
                }

            // Sets descriptors
            for (int i=0; i<MaxPath; i++) {
                if(i==7){
                    if (GGIval[i]>0)
                        this.setGGI8(GGIval[i]);
                    else
                        this.setGGI8(0);
                    return;
                }
            }

        } catch (Throwable e) {
            this.setGGI8(MISSING_VALUE);
        }
    }

    private void CalculateEA(InsilicoMolecule Mol){
        this.setEEig10ri(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setEEig10ri(MISSING_VALUE);
            return;
        }

        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            this.setEEig10ri(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = Mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setEEig10ri(MISSING_VALUE);
            return;
        }

        Matrix DataMatrix = null;

        double[][] EdgeResMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeResMat.length; i++)
            for (int j=0; j<EdgeResMat[0].length; j++)
                EdgeResMat[i][j] = EdgeAdjMat[i][j][1];

        for (int i=0; i<m.getBondCount(); i++) {
            EdgeResMat[i][i] = GetResonanceIntegral(m.getBond(i));
            if (EdgeResMat[i][i] == 0)
                EdgeResMat[i][i] = 1;

        }

        DataMatrix = new Matrix(EdgeResMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        for (int i=1; i<=15; i++) {
            if(i==10){
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    this.setEEig10ri(eigenvalues[idx]);
                else
                    this.setEEig10ri(0);
                return;
            }

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

    private void CalculateTopological(InsilicoMolecule Mol){
        this.setGNar(0);

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setGNar(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setGNar(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int[] VD = VertexDegree.getWeights(curMol, true);

            // Wiener index
            double W = 0;
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    W += TopoMat[i][j];
                }



            // Narumi indices
            double SNar=1, GNar=1;
            for (int i=0; i<nSK; i++) {
                SNar *= VD[i];
            }
            GNar = Math.pow(SNar, 1.00/nSK);
            this.setGNar(GNar);


        } catch (Throwable e) {
            log.warn(e.getMessage());
            this.setGNar(MISSING_VALUE);
        }




    }

    private void CalculateConstitutional(InsilicoMolecule Mol){
        this.setNN(0);

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setNN(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();

            int nN=0;

            for (int i=0; i<nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                if (CurAt.getSymbol().equalsIgnoreCase("N"))
                    nN++;
            }
            this.setNN(nN);
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setNN(MISSING_VALUE);
        }
    }

    private void CalculateACF(InsilicoMolecule Mol){
        DescriptorBlock acf = new ACF();
        try {
            acf.Calculate(Mol);
            this.setO_058(acf.GetByName("O-058").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setO_058(MISSING_VALUE);
        }
    }

    private void CalculateFG(InsilicoMolecule Mol){
        DescriptorBlock fg = new FG();
        try {
            fg.Calculate(Mol);
            this.setNCconj(fg.GetByName("nCconj").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setNCconj(MISSING_VALUE);
        }
    }

    private double GetResonanceIntegral(IBond bnd) {

        IAtom atA =  bnd.getAtom(0);
        IAtom atB =  bnd.getAtom(1);
        String A = atA.getSymbol();
        String B = atB.getSymbol();

        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("C") == 0)) )
            return 1.00;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("B") == 0)) ||
                ((A.compareToIgnoreCase("B") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("N") == 0)) ||
                ((A.compareToIgnoreCase("N") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.9;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("O") == 0)) ||
                ((A.compareToIgnoreCase("O") == 0) && (B.compareToIgnoreCase("C") == 0))) {
            if (bnd.getOrder() == IBond.Order.SINGLE)
                return 0.8;
            if (bnd.getOrder() == IBond.Order.DOUBLE)
                return 1.2;
        }
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("S") == 0)) ||
                ((A.compareToIgnoreCase("S") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("F") == 0)) ||
                ((A.compareToIgnoreCase("F") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Cl") == 0)) ||
                ((A.compareToIgnoreCase("Cl") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.4;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Br") == 0)) ||
                ((A.compareToIgnoreCase("Br") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.3;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("I") == 0)) ||
                ((A.compareToIgnoreCase("I") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.1;

        return 0.00;
    }
}
