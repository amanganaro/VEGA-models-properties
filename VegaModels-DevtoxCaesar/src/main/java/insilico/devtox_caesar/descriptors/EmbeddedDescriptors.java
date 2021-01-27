package insilico.devtox_caesar.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.old.weight.EState;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.devtox_caesar.descriptors.weights.BurdenMatrix;
import insilico.devtox_caesar.descriptors.weights.Mass;
import insilico.devtox_caesar.descriptors.weights.Polarizability;
import insilico.devtox_caesar.descriptors.weights.VanDerWaals;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private int MISSING_VALUE = -999;
    private double BEH1m = MISSING_VALUE;
    private double BEL3p = MISSING_VALUE;
    private double BEL1v = MISSING_VALUE;
    private double BEL8v = MISSING_VALUE;
    private double GATS1p = MISSING_VALUE;
    private double GATS2m = MISSING_VALUE;
    private double GATS3v = MISSING_VALUE;
    private double MATS1p = MISSING_VALUE;
    private double MATS4p = MISSING_VALUE;
    private double MATS4v = MISSING_VALUE;
    private double SdssC = MISSING_VALUE;
    private double SHssNH = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        try {
            CalculateAllDescriptors(Mol);
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalcuateBEH(Mol);
        CalculateAutocorrelation(Mol);
        CalculateEStates(Mol);
        CalculateAutoCorrelationPweighetd(Mol);
    }

    private void CalculateAutoCorrelationPweighetd(InsilicoMolecule Mol){
        GATS1p = 0; MATS1p = 0; MATS4p = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            GATS1p = MISSING_VALUE; MATS1p = MISSING_VALUE; MATS4p = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            GATS1p = MISSING_VALUE; MATS1p = MISSING_VALUE; MATS4p = MISSING_VALUE;
            return;
        }


        int nSK = m.getAtomCount();
        double[] w = Polarizability.getWeights(m);

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double MoranAC=0, GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 1) {
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


        this.setMATS1p(MoranAC);
        this.setGATS1p(GearyAC);

        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        MoranAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 4) {
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

        this.setMATS4p(MoranAC);


    }

    private void CalculateEStates(InsilicoMolecule Mol){
        SdssC = MISSING_VALUE;
        SHssNH = MISSING_VALUE;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SdssC = 0;
            SHssNH = 0;
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            SdssC = MISSING_VALUE;
            SHssNH = MISSING_VALUE;
            return;
        }

        // Calculation
        double SsCl=0, SsF=0, SsBr=0, SsI=0;
        double StN=0, SsssN=0, SsSH=0;
        double SdssC=0, SHssNH=0, SdssNp=0;
        double StsC=0, SsOm=0;

        double Hmax= Descriptor.MISSING_VALUE, Hmin= Descriptor.MISSING_VALUE;
        double Gmax= Descriptor.MISSING_VALUE, Gmin= Descriptor.MISSING_VALUE;

        double Ss = 0, Ms = 0;

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn("unable to get H count");
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


            // C Groups
            if (curAt.getSymbol().equalsIgnoreCase("C")) {

                if ((nBnd == 3) && (nDbl == 1) && (nSng == 2))
                    SdssC += es.getEState()[at];
            }

            // N Groups
            if (curAt.getSymbol().equalsIgnoreCase("N")) {
                if ((nBnd == 2) && (nSng == 2) && (nH == 1))
                    SHssNH += es.getHEState()[at];
            }

        }


        this.setSdssC(SdssC);
        this.setSHssNH(SHssNH);
    }

    private void CalculateAutocorrelation(InsilicoMolecule Mol) {
        GATS2m = 0; GATS3v = 0; MATS4v = 0;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            GATS2m = MISSING_VALUE; GATS3v = MISSING_VALUE;MATS4v = MISSING_VALUE;
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            GATS2m = MISSING_VALUE; GATS3v = MISSING_VALUE;MATS4v = MISSING_VALUE;
            return;
        }


        int nSK = m.getAtomCount();

        // m
        double[] w = Mass.getWeights(m);
        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        // Calculates autocorrelations

        double GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 2) {
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


        this.setGATS2m(GearyAC);

        // v
        w = VanDerWaals.getWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        GearyAC=0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 3) {
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

        this.setGATS3v(GearyAC);

        double MoranAC =0;
        denom = 0; delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 4) {

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


        this.setMATS4v(MoranAC);
    }

    private void CalcuateBEH(InsilicoMolecule Mol){
        BEH1m = 0; BEL1v = 0; BEL3p = 0; BEL8v = 0;

        InsilicoMolecule HMol;
        IAtomContainer m;
        try {
            HMol = (InsilicoMolecule) Mol.Clone();
            HMol.SetExplicitHydrogen(true);
            m = HMol.GetStructure();
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            BEH1m = MISSING_VALUE; BEL1v = MISSING_VALUE; BEL3p = MISSING_VALUE; BEL8v = MISSING_VALUE;
            return;
        }

        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = BurdenMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            BEH1m = MISSING_VALUE; BEL1v = MISSING_VALUE; BEL3p = MISSING_VALUE; BEL8v = MISSING_VALUE;
            return;
        }

        int nSK = m.getAtomCount();

        // m
        double[] w = Mass.getWeights(m);
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
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valH = 0;
            } else {
                valH = eigenvalues[eigenvalues.length-1-i];
            }
            if( i == 0)
                this.setBEH1m(valH);
        }

        // p
        w = Polarizability.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valL = 0;
            } else {
                valL = eigenvalues[i];
            }
            if( i == 2)
                this.setBEL3p(valL);
        }

        // v
        w = VanDerWaals.getWeights(m);
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        DataMatrix = new Matrix(BurdenMat);
        ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);

        for (int i=0; i<8; i++) {
            double valH, valL;
            if (i>(eigenvalues.length-1)) {
                valL = 0;
            } else {
                valL = eigenvalues[i];
            }
            if( i == 0)
                this.setBEL1v(valL);
            if( i == 7)
                this.setBEL8v(valL);
        }

    }






}
