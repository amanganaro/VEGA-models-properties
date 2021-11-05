package insilico.logk.descriptors;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.descriptor.blocks.*;
import insilico.descriptor.blocks.logP.ALogP;
import insilico.descriptor.blocks.logP.MLogP;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtom;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;

@Slf4j
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    public double ALogP = MISSING_VALUE;
    public double P_VSA_i_2 = MISSING_VALUE;
    public double MLogP = MISSING_VALUE;
    public double P_VSA_p_3 = MISSING_VALUE;
    public double Eta_betaP = MISSING_VALUE;
    public double CATS2D_00_LL = MISSING_VALUE;
    public double C = MISSING_VALUE;
    public double PCD = MISSING_VALUE;
    public double Ui = MISSING_VALUE;
    public double N = MISSING_VALUE;
    public double totalcharge = MISSING_VALUE;
    public double CATS2D_00_PP = MISSING_VALUE;
    public double C_024 = MISSING_VALUE;
    public double nCsp2 = MISSING_VALUE;

    public double MLOGP2 = MISSING_VALUE;
    public double GATS1i = MISSING_VALUE;
    public double SpMax2_Bh_P = MISSING_VALUE;
    public double nBm = MISSING_VALUE;
    public double MATS5e = MISSING_VALUE;
    public double AMW = MISSING_VALUE;
    public double F01_C_N = MISSING_VALUE;
    public double T_O_O = MISSING_VALUE;
    public double J_D_DT = MISSING_VALUE;
    public double SpMax_AEA_dm = MISSING_VALUE;

    public double ExperimentalValue;

    public double[] getDescriptors(){
        return new double[] {
                ALogP, P_VSA_i_2, MLogP, P_VSA_p_3, Eta_betaP, CATS2D_00_LL, C,
                PCD, Ui, N, totalcharge, CATS2D_00_PP, C_024, nCsp2, MLOGP2, GATS1i, SpMax2_Bh_P,
                nBm, MATS5e, AMW, F01_C_N, T_O_O, J_D_DT, SpMax_AEA_dm};
    }


    public EmbeddedDescriptors(InsilicoMolecule mol, boolean fromFile) throws MalformedURLException {
        if(fromFile)
            SearchDescriptors(mol);
        else CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) {
        CalculateLogP(mol);
        CalculateP_VSA(mol);
        CalculateEtaBeta(mol);
        CalculateCATS2D(mol);
        CalculateConstitutional(mol);
        CalculateWAP(mol);
        CalculateACF(mol);
        CalculateAutoCorrelation(mol);
        CalculateSP(mol);
        CalculateSPMax(mol);
        CalculateAtomPairs(mol);
        CalculateFG(mol);
        CalculateJDDT(mol);
        CalculateNcs(mol);
    }

    private void CalculateNcs(InsilicoMolecule mol) {
        DescriptorBlock block = new FunctionalGroups();
        try {
            block.Calculate(mol);
            nCsp2 = block.GetByName("nCs").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateJDDT(InsilicoMolecule mol) {
        DescriptorBlock block = new Matrices2D();
        try {
            block.Calculate(mol);
            J_D_DT = block.GetByName("SpAbs_D/Dt").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }


    private void CalculateFG(InsilicoMolecule mol) {

//        DescriptorBlock block = new FunctionalGroups();
//        try {
//            block.Calculate(mol);
//            nCsp2 = block.GetByName("nCs").getValue();
//        } catch (DescriptorNotFoundException ex){
//            log.warn(ex.getMessage());
//        }

    }

    private void CalculateAtomPairs(InsilicoMolecule mol) {

        DescriptorBlock block = new AtomPairs2D();
        try {
            block.Calculate(mol);
            F01_C_N = block.GetByName("F01[C-N]").getValue();
            T_O_O = block.GetByName("T(O..O)").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }


    }

    private void CalculateSPMax(InsilicoMolecule mol) {

        DescriptorBlock block = new EdgeAdjacency();
        try {
            block.Calculate(mol);
            SpMax_AEA_dm = block.GetByName("SpMaxA_EA(dm)").getValue();
        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


    }

    private void CalculateSP(InsilicoMolecule mol) {
        BurdenEigenvalues block = new BurdenEigenvalues();

        try {
            block.Calculate(mol);
            SpMax2_Bh_P = block.GetByName("SpMax2_Bh(p)").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateAutoCorrelation(InsilicoMolecule mol) {

        AutoCorrelation block = new AutoCorrelation();
        try {
            block.Calculate(mol);
            GATS1i = block.GetByName("GATS1i").getValue();
            MATS5e = block.GetByName("MATS5e").getValue();

        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateACF(InsilicoMolecule mol) {
        AtomCenteredFragments block = new AtomCenteredFragments();

        try {
            block.Calculate(mol);
            C_024 = block.GetByName("C-024").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateWAP(InsilicoMolecule mol) {
        DescriptorBlock block = new WalkAndPath();

        try {
            block.Calculate(mol);
            PCD = block.GetByName("PCD").getValue();


//            MoleculePaths paths = new MoleculePaths(mol);
//            double TPC = paths.Total_Path_Count;
//            double piID = paths.IDpi;
//            TPC = Math.log(1 + TPC);
//            piID = Math.log(1 + piID);

//            PCD =  Math.round(piID - TPC);
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateConstitutional(InsilicoMolecule mol) {

        DescriptorBlock block = new Constitutional();
        block.Calculate(mol);
        try {

            N = block.GetByName("NPerc").getValue();
            C = block.GetByName("CPerc").getValue();
            nBm = block.GetByName("nBm").getValue();
            AMW = block.GetByName("AMW").getValue();


            double scbo = block.GetByName("SCBO").getValue();
            double nBo = block.GetByName("nBo").getValue();

            Ui = Math.log(1 + scbo - nBo) / Math.log(2);

            totalcharge = 0;
            for (IAtom at : mol.GetStructure().atoms())
                totalcharge += at.getFormalCharge();

        } catch (DescriptorNotFoundException | InvalidMoleculeException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateCATS2D(InsilicoMolecule mol) {

        DescriptorBlock block = new Cats2D();

        try {
            block.Calculate(mol);

            CATS2D_00_PP = block.GetByName("CATS2D_00_PP").getValue();
            CATS2D_00_LL = block.GetByName("CATS2D_00_LL").getValue();



        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateEtaBeta(InsilicoMolecule mol) {
        DescriptorBlock block = new EtaIndices();

        try {
            block.Calculate(mol);
            Eta_betaP = block.GetByName("Eta_betaP").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }

    }

    private void CalculateP_VSA(InsilicoMolecule mol) {
        DescriptorBlock block = new P_VSA();

        try {
            block.Calculate(mol);
            P_VSA_p_3 = block.GetByName("P_VSA_p_3").getValue();
            P_VSA_i_2 = block.GetByName("P_VSA_i_2").getValue();
        } catch (DescriptorNotFoundException ex){
            log.warn(ex.getMessage());
        }
    }

    private void CalculateLogP(InsilicoMolecule mol) {
        DescriptorBlock block = new MLogP();
        try {
            block.Calculate(mol);
            MLogP = block.GetByName("MLogP").getValue();
        } catch (DescriptorNotFoundException ex) {
            log.warn(ex.getMessage());
        }

        block = new ALogP();
        try {
            block.Calculate(mol);
            ALogP = block.GetByName("ALogP").getValue();
        } catch (DescriptorNotFoundException ex) {
            log.warn(ex.getMessage());
        }

        MLOGP2 = Math.pow(MLogP, 2);
    }

    private void SearchDescriptors(InsilicoMolecule mol) throws MalformedURLException {

        URL url = new URL("file:///" + System.getProperty("user.dir") + "/VegaModels-LogK/src/main/resources/data/dataset_logk.csv");

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))) {

            // skip the first row (header)
            br.readLine();

            String line;
            while ((line = br.readLine()) != null){
                String[] lineArray = line.split("\t");

                // compare
                if(mol.GetSMILES().equals(SmilesMolecule.Convert(lineArray[2].trim()).GetSMILES())){

                    ExperimentalValue = Double.parseDouble(lineArray[3]);

                    ALogP = Double.parseDouble(lineArray[7]);
                    P_VSA_i_2 = Double.parseDouble(lineArray[8]);
                    MLogP = Double.parseDouble(lineArray[9]);
                    P_VSA_p_3 = Double.parseDouble(lineArray[10]);
                    Eta_betaP = Double.parseDouble(lineArray[11]);
                    CATS2D_00_LL = Double.parseDouble(lineArray[12]);
                    C = Double.parseDouble(lineArray[13]);
                    PCD = Double.parseDouble(lineArray[14]);
                    Ui = Double.parseDouble(lineArray[15]);
                    N = Double.parseDouble(lineArray[16]);
                    totalcharge = Double.parseDouble(lineArray[17]);
                    CATS2D_00_PP = Double.parseDouble(lineArray[18]);
                    C_024 = Double.parseDouble(lineArray[19]);
                    nCsp2 = Double.parseDouble(lineArray[20]);
                    MLOGP2 = Double.parseDouble(lineArray[21]);
                    GATS1i = Double.parseDouble(lineArray[22]);
                    SpMax2_Bh_P = Double.parseDouble(lineArray[23]);
                    nBm = Double.parseDouble(lineArray[24]);
                    MATS5e = Double.parseDouble(lineArray[25]);
                    AMW = Double.parseDouble(lineArray[26]);
                    F01_C_N = Double.parseDouble(lineArray[27]);
                    T_O_O = Double.parseDouble(lineArray[28]);
                    J_D_DT = Double.parseDouble(lineArray[29]);
                    SpMax_AEA_dm = Double.parseDouble(lineArray[30]);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
            log.warn(e.getMessage());
        }

    }


}
