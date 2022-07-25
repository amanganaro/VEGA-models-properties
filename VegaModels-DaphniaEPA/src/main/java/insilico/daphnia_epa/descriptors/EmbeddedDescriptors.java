package insilico.daphnia_epa.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.daphnia_epa.descriptors.weights.*;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;

@Log4j
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;
    public double StN = MISSING_VALUE;
    public double SsSH = MISSING_VALUE;
    public double nROH = MISSING_VALUE;
    public double nArOH = MISSING_VALUE;
    public double HMax = MISSING_VALUE;
    public double MDE_N_33 = MISSING_VALUE;
    public double BEH1m = MISSING_VALUE;
    public double BEH1p = MISSING_VALUE;
    public double Mv = MISSING_VALUE;
    public double Mw = MISSING_VALUE;
    public double MATS1m = MISSING_VALUE;
    public double MATS1e = MISSING_VALUE;
    public double GATS3m = MISSING_VALUE;
    public double AMR = MISSING_VALUE;


    public EmbeddedDescriptors(InsilicoMolecule Mol) {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateEStates(Mol);
        CalculateFG(Mol);
        CalculateDE(Mol);
        CalculateBEH(Mol);
        CalculateConstitutional(Mol);
        CalculateAutocorrelation(Mol);
        CalculateALogP(Mol);
    }

    private void CalculateALogP(InsilicoMolecule Mol){

        AMR = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            AMR = MISSING_VALUE;
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF = new AtomCenteredFragments();
        CurACF.Calculate(Mol);


        double[] Frags = CurACF.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            AMR += Frags[i] * GhoseCrippenWeights.GetMolarRefractivity(i);
        }
        
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol){
        MATS1e = 0; MATS1m = 0; GATS3m = 0;
        int lag1 = 1;
        int lag3 = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            MATS1e = MISSING_VALUE; MATS1m = MISSING_VALUE; GATS3m = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MATS1e = MISSING_VALUE; MATS1m = MISSING_VALUE; GATS3m = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w = Electronegativity.getWeights(m);
        // Calculates weights averages
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag1) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MATS1e = 1;
            } else {
                MATS1e = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        w = Mass.getWeights(m);
        // Calculates weights averages
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        MoranAC=0;
        denom = 0;
        delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag1) {
                    MoranAC += (w[i] - wA) * (w[j] - wA);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                MATS1m = 1;
            } else {
                MATS1m = ((1 / delta) * MoranAC) / ((1 / ((double)nSK)) * denom);
            }
        }

        w = Mass.getWeights(m);
        // Calculates weights averages
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag3) {

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

        GATS3m = GearyAC;

    }

    private void CalculateConstitutional(InsilicoMolecule Mol){
        DescriptorBlock block = new Constitutional();
        block.Calculate(Mol);
        try {
            Mw = block.GetByName("MW").getValue();
            Mv = block.GetByName("Mv").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateBEH(InsilicoMolecule Mol){
        BEH1m = 0; BEH1p = 0;

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            BEH1m = MISSING_VALUE; BEH1p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat_m;
        double[][] BurdenMat_p;
        try {
            BurdenMat_m = BurdenMatrix.getMatrix(m);
            BurdenMat_p = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            BEH1m = MISSING_VALUE; BEH1p = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        double[] w_p = Polarizability.getWeights(m);
        double[] w_m = Mass.getWeights(m);

        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++) {
            if (w_p[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
            if (w_m[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
        }

        for (int i=0; i<nSK; i++) {
            BurdenMat_m[i][i] = w_m[i];
            BurdenMat_p[i][i] = w_p[i];
        }

        Matrix DataMatrixm = new Matrix(BurdenMat_m);
        Matrix DataMatrixp = new Matrix(BurdenMat_p);
        double[] eigenvaluesm;
        double[] eigenvaluesp;
        EigenvalueDecomposition edm = new EigenvalueDecomposition(DataMatrixm);
        EigenvalueDecomposition edp = new EigenvalueDecomposition(DataMatrixp);
        eigenvaluesm = edm.getRealEigenvalues();
        eigenvaluesp = edp.getRealEigenvalues();
        Arrays.sort(eigenvaluesm);
        Arrays.sort(eigenvaluesp);

        for (int i=0; i<10; i++) {
            double valH;
            if (i==0){
                if (i>(eigenvaluesm.length-1)) {
                    BEH1m = 0;
                } else {
                    BEH1m = eigenvaluesm[eigenvaluesm.length-1-i];
                }

                if (i>(eigenvaluesp.length-1)) {
                    BEH1p = 0;
                } else {
                    BEH1p = eigenvaluesp[eigenvaluesp.length-1-i];
                }
            }
        }

    }


    private void CalculateDE(InsilicoMolecule Mol){
        MDE_N_33 = 0;
        int vd03 = 3;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            MDE_N_33 = MISSING_VALUE;
            return;
        }

        // Get matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            MDE_N_33 = MISSING_VALUE;
            return;
        }

        int[] VD = VertexDegree.getWeights(m, true);
        int nSK = m.getAtomCount();

        double n = 0;
        double d = 1;

        for (int i=0; i<nSK; i++)
            if ("N".equalsIgnoreCase(m.getAtom(i).getSymbol()))
                for (int j=i+1; j<nSK; j++)
                    if ("N".equalsIgnoreCase(m.getAtom(j).getSymbol()))
                        if  ( ((VD[i]==vd03) && (VD[j]==vd03)) || ((VD[i]==vd03) && (VD[j]==vd03)) ) {
                            n++;
                            d *= TopoMatrix[i][j];
                        }

//        double res = 0;
        if (n > 0) {
            d = Math.pow(d,1.0/(2.0*n));
            MDE_N_33 = n/(d*d);
        }

    }


    private void CalculateFG(InsilicoMolecule Mol){
        DescriptorBlock descriptorBlock = new FunctionalGroups();
        descriptorBlock.Calculate(Mol);
        try {
            nArOH = descriptorBlock.GetByName("nArOH").getValue();
            nROH = descriptorBlock.GetByName("nROH").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateEStates(InsilicoMolecule Mol){
        StN = 0; SsSH = 0; HMax = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            StN = MISSING_VALUE; SsSH = MISSING_VALUE; HMax = MISSING_VALUE;
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            StN = MISSING_VALUE; SsSH = MISSING_VALUE; HMax = MISSING_VALUE;
            return;
        }

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("unable to get H count");
            }

            // formal charge
            int Charge;
            try {
                Charge = curAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }

            // Count bonds
            int nBnd=0, nSng = 0, nDbl = 0, nTri = 0, nAr=0;
            for (IBond b : m.getConnectedBondsList(curAt)) {
                if (b.getFlag(CDKConstants.ISAROMATIC)) {
                    nAr++;
                    nBnd++;
                    continue;
                }
                if (b.getOrder() == IBond.Order.SINGLE) {
                    nSng++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.DOUBLE) {
                    nDbl++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nTri++;
                    nBnd++;
                }
            }

            if (nH>0) {
                HMax = (HMax == MISSING_VALUE) ? es.getHEState()[at] : (Math.max(es.getHEState()[at], HMax));
            }

            // S groups
            if (curAt.getSymbol().equalsIgnoreCase("S")) {

                if ((nBnd == 1) && (nSng == 1) && (nH == 1))
                    SsSH += es.getEState()[at];
            }

            // N groups
            if (curAt.getSymbol().equalsIgnoreCase("N")) {

                if ((nBnd == 1) && (nTri == 1) )
                    StN += es.getEState()[at];
            }

        }

    }
}
