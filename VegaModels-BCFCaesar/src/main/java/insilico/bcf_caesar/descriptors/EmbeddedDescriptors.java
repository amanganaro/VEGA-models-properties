package insilico.bcf_caesar.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.bcf_caesar.descriptors.weights.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;


public class EmbeddedDescriptors {
    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;
    public double MLogP = MISSING_VALUE;
    public double X0sol = MISSING_VALUE;
    public double MATS5v = MISSING_VALUE;
    public double GATS5v = MISSING_VALUE;
    public double BEH2p = MISSING_VALUE;
    public double AEige = MISSING_VALUE;
    public double Cl_089 = MISSING_VALUE;
    public double ssCl = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
        MLogP = descriptorMLogP.MLogP;

        CalculateX0sol(Mol);
        CalculateMats(Mol);
        CalculateBEH2p(Mol);
        CalculateAEige(Mol);
        CalculateCl089(Mol);
        CalculatessCL(Mol);
    }

    private void CalculatessCL(InsilicoMolecule Mol){

        ssCl = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            ssCl = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            ssCl = MISSING_VALUE;
            return;
        }

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Halo atoms
            if (curAt.getSymbol().equalsIgnoreCase("Cl")) {
                ssCl += es.getEState()[at];
            }

        }
    }

    private void CalculateCl089(InsilicoMolecule Mol) throws DescriptorNotFoundException {
        ACF acf = new ACF();
        acf.Calculate(Mol);
        Cl_089 = acf.GetByName("Cl-089").getValue();
    }

    private void CalculateAEige(InsilicoMolecule Mol){
        AEige = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            AEige = MISSING_VALUE;
            return;
        }

        double[][] ConnMatrix;
        try {
            ConnMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            AEige = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w = Electronegativity.getWeights(m);;
        double refW = Electronegativity.GetElectronegativity("C");

        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++)
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;

        // Builds matrix
        double[][] EigMat = new double[nSK][nSK];
        for (int i=0; i<nSK; i++)
            for (int j=0; j<nSK; j++) {

                if (i==j) {

                    EigMat[i][j] = 1 - (refW / w[i]);

                } else {

                    // builds shortest path between i and j
                    IAtom at1 = m.getAtom(i);
                    IAtom at2 =  m.getAtom(j);
                    ShortestPaths sp = new ShortestPaths(m, at1);
                    List<IAtom> Path = Arrays.asList(sp.atomsTo(at2));

                    double val = 0;
                    for (int k=0; k<(Path.size()-1); k++) {
                        int a1 = m.indexOf(Path.get(k));
                        int a2 = m.indexOf(Path.get(k + 1));
                        double bond = ConnMatrix[a1][a2];
                        val += (1 / bond) * (Math.pow(refW, 2) / (w[a1] * w[a2]) );
                    }

                    EigMat[i][j] = val;

                }
            }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(EigMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<eigenvalues.length; i++) {
            AEige += Math.abs(eigenvalues[i]);
        }


    }

    private void CalculateBEH2p(InsilicoMolecule Mol) {

        BEH2p = 0;


        int MaxEig = 8;

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            log.warn(e.getMessage());
            BEH2p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat_p;
        BurdenMat_p = BurdenMatrix.getMatrix(m);

        int nSK = m.getAtomCount();
        double[] w_p = Polarizability.getWeights(m);

        for (int i=0; i<nSK; i++) {
            BurdenMat_p[i][i] = w_p[i];
        }

        Matrix DataMatrix_p = new Matrix(BurdenMat_p);
        double[] eigenvalues_p;
        EigenvalueDecomposition ed_p = new EigenvalueDecomposition(DataMatrix_p);
        eigenvalues_p = ed_p.getRealEigenvalues();
        Arrays.sort(eigenvalues_p);

        // p
        for (int i=0; i<MaxEig; i++) {
            if(i == 1) {
                double valH;
                if (i>(eigenvalues_p.length-1)) {
                    valH = 0;
                } else {
                    valH = eigenvalues_p[eigenvalues_p.length-1-i];
                }
                BEH2p = valH;
            }

        }

    }

    private void CalculateMats(InsilicoMolecule Mol){
        MATS5v = 0;
        GATS5v = 0;

        int PARAMETER_LAG= 5;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            MATS5v = MISSING_VALUE;
            GATS5v = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS5v = MISSING_VALUE;
            GATS5v = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double w[] = VanDerWaals.getWeights(m);
        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++) {
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0, GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == PARAMETER_LAG) {
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

        GATS5v = GearyAC;
        MATS5v = MoranAC;
    }

    private void CalculateX0sol(InsilicoMolecule Mol) {
        X0sol = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(e.getMessage());
            X0sol = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] ConnAugMatrix;
        try {
            ConnAugMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            X0sol = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        int[] VD = VertexDegree.getWeights(m, true);
        double[] ValenceVD = ValenceVertexDegree.getWeights(m);
        int[] Qnumbers = QuantumNumber.getWeights(m);


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


        for (int i=0; i<nSK; i++) {
            if (ConnAugMatrix[i][i] == 9)
                continue; // F not taken into account
            // path 0
            X0sol += 0.5 * Qnumbers[i] * Math.pow(VD[i], -0.5);


        }
    }


}
