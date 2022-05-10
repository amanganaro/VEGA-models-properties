package insilico.ppara_up.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;

import insilico.ppara_up.utils.MoleculePaths;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;

@Slf4j
public class EmbeddedDescriptors {

    private final double MISSING_VALUE = -999;

    public double TI2_L = MISSING_VALUE;
    public double PW4 = MISSING_VALUE;
    public double CATS2D_07_NL = MISSING_VALUE;
    public double ExperimentalValue = MISSING_VALUE;

    public double[] getDescriptors(){
        return new double[]{TI2_L, PW4, CATS2D_07_NL};
    }

    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        if(fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
        System.out.println();
    }


    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateTI2(mol);
        CalculatePW4(mol);
        CalculateCATS2D(mol);
    }

    private void CalculateTI2(InsilicoMolecule mol) {
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure for: " + mol.GetSMILES());
            return;
        }

        try {
            // Adj matrix available for calculations
            int[][] AdjMatrix = mol.GetMatrixAdjacency();

            int nSK = curMol.getAtomCount();
            int nBo = curMol.getBondCount();
            double[][] Mat = new double[nSK][nSK];

            int[][] LapMatrix = mol.GetMatrixLaplace();
            for (int i=0; i<nSK; i++)
                for (int j=0; j<nSK; j++)
                    Mat[i][j] = LapMatrix[i][j];

            Matrix DataMatrix = new Matrix(Mat);

            double[] eigenvalues = new double[0];

            try {
                EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
                eigenvalues = ed.getRealEigenvalues();
                Arrays.sort(eigenvalues);
            } catch (Throwable e) {
                log.warn("Unable to calculate eigenvalue: " + e.getMessage());
            }

            double MIN_EIG_VALUE = Math.pow(1.5, -45);

            double LastEigen = 0;
            int matrix_dim = eigenvalues.length;

            for (int i=eigenvalues.length-1; i>=0; i--) {
                double val = eigenvalues[i];
//                for (double val : eigenvalues) {
                if (Math.abs(val) < MIN_EIG_VALUE)
                    continue;
                if (val > 0) {
                    LastEigen = val;
                }
            }
            if (LastEigen > 0)
                TI2_L = 4.0 / (matrix_dim * LastEigen);
            else
                TI2_L = Descriptor.MISSING_VALUE;

        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }




    }

    private void CalculatePW4(InsilicoMolecule mol) {
        try {

            MoleculePaths paths = new MoleculePaths(mol);
            PW4 = paths.Pws[3];

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateCATS2D(InsilicoMolecule mol) {

        String[] TYPE_N = { "N", ""};
        String[] TYPE_L = { "L", ""};

        String[][][] AtomCouples = {
                {TYPE_N, TYPE_L}
        };

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure for: " + mol.GetSMILES());
            return;
        }

        int nSK = curMol.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        double[][] ConnAugMatrix = null;
        try {
            ConnAugMatrix = mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return;
        }

        ArrayList<String>[] CatsTypes = setCatsAtomType(curMol, ConnAugMatrix);


        for (String[][] atomCouple : AtomCouples) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i = 0; i < nSK; i++) {
                if (isIn(CatsTypes[i], atomCouple[0][0])) {
                    for (int j = i; j < nSK; j++) {
//                        if (i==j) continue;
                        if (isIn(CatsTypes[j], atomCouple[1][0])) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i = 0; i < desc.length; i++)
                if(i==7)
                    CATS2D_07_NL = desc[i];

        }
    }


    private void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {
        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-PPARa_up/src/main/resources/data/dataset.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            br.readLine();
            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[2].trim()).GetSMILES())) {


                    ExperimentalValue = Double.parseDouble(lineArray[5]);

                    //6 7 8
                    TI2_L = Double.parseDouble(lineArray[7]);
                    PW4 = Double.parseDouble(lineArray[8]);
                    CATS2D_07_NL = Double.parseDouble(lineArray[9]);
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }

    }

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }

    private ArrayList<String>[] setCatsAtomType(IAtomContainer m, double[][]ConnAugMatrix) {

        int nSK = m.getAtomCount();
        ArrayList[] AtomTypes = new ArrayList[nSK];

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
            } catch (Exception e) {
                log.warn("Unable to count H");
            }

            // counters
            int nSglBnd = 0, nOtherBnd = 0, VD = 0;
            int nC = 0, nDblO = 0, nOtherNonOBond=0, nSglOH = 0;
            for (int j=0; j<nSK; j++) {
                if (j==i) continue;
                if (ConnAugMatrix[i][j]>0) {

                    VD++;

                    if (ConnAugMatrix[j][j] == 6)
                        nC++;

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
                        nOtherBnd++;
                        if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
                            nDblO++;
                        else
                            nOtherNonOBond++;
                    }
                }

            }


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

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (H == 1) || (H ==2) )
                    tD = true;

            }

            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                    (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                    (CurAt.getFormalCharge() == 0) )  {

                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherNonOBond == 0) )
                    tN = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I"))
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                if ( VD == nC)
                    tL = true;
            }

            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                if ( (nC == 2) && ( (VD+H) == 2 ))
                    tL = true;
            }


            // Sets final types
            if (tN) AtomTypes[i].add("N");
            if (tL) AtomTypes[i].add("L");

        }

        return AtomTypes;
    }


}
