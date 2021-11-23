package insilico.ahr_up.descriptor;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.iBasicWeight;
import insilico.core.descriptor.blocks.weights.other.WeightsIState;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.descriptor.blocks.AutoCorrelation;
import insilico.descriptor.blocks.BurdenEigenvalues;
import insilico.descriptor.blocks.EdgeAdjacency;
import insilico.descriptor.blocks.Matrices2D;
import insilico.descriptor.localization.StringSelectorDescriptors;
import javassist.runtime.Desc;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;

@Slf4j
public class EmbeddedDescriptors {

    private final double MISSING_VALUE = -999;

    public double ATSC5m = MISSING_VALUE;
    public double SpMin3_Bh_v = MISSING_VALUE;
    public double SpMaxA_B_s = MISSING_VALUE;
    public double ChiA_D = MISSING_VALUE;
    public double ATS7s = MISSING_VALUE;
    public double ATSC7e = MISSING_VALUE;
    public double Si = MISSING_VALUE;

    public double[] getDescriptors(){
        return new double[]{
                ATSC5m,  SpMin3_Bh_v, SpMaxA_B_s, ChiA_D, ATS7s, ATSC7e, Si
        };
    }

    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        if(fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateATS(mol);
        CalculateSpMax(mol);
        CalculateSpMin(mol);
        CalculateChiA(mol);
        CalculateSi(mol);
    }

    private void CalculateSi(InsilicoMolecule mol) {
        IAtomContainer curMol = null;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn(String.format(StringSelectorCore.getString("descriptors_invalid_structure"), mol.GetSMILES()));
        }

        try {
            int nSK = curMol.getAtomCount();
            int[] H = new int[nSK];


            iBasicWeight ws = new WeightsIonizationPotential();

            double[] weights = ws.getScaledWeights(curMol);
            double weightH = ws.getScaledWeight("H");

            double sum = 0;
            for (int i=0; i<nSK; i++) {
                if (weights[i] == Descriptor.MISSING_VALUE) {
                    sum = Descriptor.MISSING_VALUE;
                    break;
                } else {
                    // all values INCLUDING MW are scaled on carbon
                    sum += weights[i];
                    if (H[i] > 0)
                        sum += weightH * H[i];
                }
            }

            Si = sum;
        } catch (Throwable e){
            log.warn(String.format(StringSelectorCore.getString("descriptors_unable_calculate"), "Constitutional Descriptors", e.getMessage()));

        }
    }

    private void CalculateChiA(InsilicoMolecule mol) {
        DescriptorBlock block = new Matrices2D();
        try {
            block.Calculate(mol);
            ChiA_D = block.GetByName("ChiA_D").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        // CHIA_D

        try {
            IAtomContainer curMol;
            try {
                curMol = mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
                return;
            }

            // Adj matrix available for calculations
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = new double[nSK][nSK];

            int[][] TopoMatrix = mol.GetMatrixTopologicalDistance();
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++)
                    Mat[i][j] = TopoMatrix[i][j];

            double[] VS = new double[nSK];
            for (int i=0; i<nSK; i++) {
                VS[i] = 0;
                for (int j=0; j<nSK; j++)
                    VS[i] += Mat[i][j];
            }

            double Chi = 0;
            double Wi = 0;
            for (int i=0; i<nSK; i++)
                Wi += Mat[i][i];

            for (int i=0; i<(nSK-1); i++)
                for (int j=(i+1); j<nSK; j++) {
                    Wi += Mat[i][j];
                    if (AdjMatrix[i][j] > 0)
                        Chi += Math.pow( (VS[i] * VS[j]) , -0.5);
                }

            ChiA_D = Chi / (double)nBo;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateSpMax(InsilicoMolecule mol) {
        DescriptorBlock block = new Matrices2D();
        try {
            block.Calculate(mol);
            SpMaxA_B_s = block.GetByName("SpMaxA_B(s)").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        try {
            IAtomContainer curMol;
            try {
                curMol = mol.GetStructure();
            } catch (InvalidMoleculeException e) {
                log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
                return;
            }

            // Adj matrix available for calculations
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();

            double[][] Mat = new double[nSK][nSK];
            Mat = mol.GetMatrixBurden();
            double[] w = (new WeightsIState()).getWeights(curMol, false);
            for (int i=0; i<nSK; i++)
                Mat[i][i] = w[i];

            Matrix DataMatrix = new Matrix(Mat);
            double[] eigenvalues;

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn(StringSelectorDescriptors.getString("unable_eigenvalue") + e.getMessage());
                return;
            }

            double EigMax = eigenvalues[0];
            double EigMin = eigenvalues[0];
            double EigExpSum = 0;
            for (double val : eigenvalues) {
                if (val > EigMax)
                    EigMax = val;
                if (val< EigMin)
                    EigMin = val;
                EigExpSum += Math.exp(val);
            }
            SpMaxA_B_s = EigMax / (double) nSK;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateSpMin(InsilicoMolecule mol) {
        DescriptorBlock block = new BurdenEigenvalues();
        try {
            block.Calculate(mol);
            SpMin3_Bh_v = block.GetByName("SpMin3_Bh(v)").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateATS(InsilicoMolecule mol) {

        DescriptorBlock block = new AutoCorrelation();
        try {
            block.Calculate(mol);
            ATS7s = block.GetByName("ATS7s").getValue();
            ATSC5m = block.GetByName("ATSC5m").getValue();
            ATSC7e = block.GetByName("ATSC7e").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

        // s
        IAtomContainer m;
        try {
            IAtomContainer orig_m = mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            log.warn(StringSelectorDescriptors.getString("invalid_structure") + mol.GetSMILES());
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            return;
        }

        int nSK = m.getAtomCount();
        double[] w = new WeightsIState().getWeights(m, true);
        for (int i=0; i<nSK; i++) {
            if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w[i] = 1;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        int lag = 7;
        double AC=0;

        for (int i=0; i<nSK; i++) {


            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag) {
                    AC += w[i] * w[j];
                }
        }

        AC /= 2.0;
        ATS7s = Math.log(1 + AC);

        // e
        w = new WeightsElectronegativity().getScaledWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        lag = 7;
        double ACS=0;

        for (int i=0; i<nSK; i++) {



            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag) {
                    ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                }
        }

        ACS /= 2.0;
        ATSC7e = ACS;


        w = new WeightsMass().getScaledWeights(m);
        wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        lag = 5;
        ACS=0;

        for (int i=0; i<nSK; i++) {



            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == lag) {
                    ACS += Math.abs((w[i]-wA) * (w[j]-wA));
                }
        }

        ACS /= 2.0;
        ATSC5m = ACS;

    }

    private void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {
        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-AhR_up/src/main/resources/data/dataset.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            br.readLine();
            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[1].trim()).GetSMILES())) {

                    //6 7 8
                    ATSC5m = Double.parseDouble(lineArray[2]);
                    SpMin3_Bh_v = Double.parseDouble(lineArray[3]);
                    SpMaxA_B_s = Double.parseDouble(lineArray[4]);
                    ChiA_D = Double.parseDouble(lineArray[5]);
                    ATS7s = Double.parseDouble(lineArray[6]);
                    ATSC7e = Double.parseDouble(lineArray[7]);
                    Si = Double.parseDouble(lineArray[8]);
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }


    }


}
