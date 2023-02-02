package insilico.vapour_pressure.descriptors;

import insilico.core.exception.GenericFailureException;
import insilico.core.molecule.InsilicoMolecule;
//import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;

    public double H_050 = MISSING_VALUE;
    public double F02_F_F = MISSING_VALUE;
    public double SM02_AEAdm = MISSING_VALUE;
    public double SpMax_AEAed = MISSING_VALUE;
    public double ATS2m = MISSING_VALUE;
    public double P_VSA_PPP_A = MISSING_VALUE;
    public double NPerc = MISSING_VALUE;
    public double piID = MISSING_VALUE;
    public double CATS2D_07_DD = MISSING_VALUE;
    public double CATS2D_01_DN = MISSING_VALUE;
    public double N_072 = MISSING_VALUE;
    public double nAmideE_ = MISSING_VALUE;
    public double X3v = MISSING_VALUE;
    public double TI1_L = MISSING_VALUE;
    public double SpMaxA_Dzm = MISSING_VALUE;
    public double BB_SA12 = MISSING_VALUE;
    public double nCONN = MISSING_VALUE;
    public double CATS2D_07_DL = MISSING_VALUE;
    public double B03_F_F = MISSING_VALUE;
    public double SaaNH = MISSING_VALUE;
    public double ED_3 = MISSING_VALUE;
    public double CATS2D_09_DA = MISSING_VALUE;
    public double LogP = MISSING_VALUE;
    public double nRNH2 = MISSING_VALUE;
    public double F07_C_S = MISSING_VALUE;
    public double F10_C_O = MISSING_VALUE;
    public double CATS2D_09_AL = MISSING_VALUE;
    public double BB_SA31c = MISSING_VALUE;
    public double SpMin5_Bhm = MISSING_VALUE;
    public double X1Av = MISSING_VALUE;
    public double N_068 = MISSING_VALUE;
    public double nRNR2 = MISSING_VALUE;
    public double F09_O_O = MISSING_VALUE;
    public double IC1 = MISSING_VALUE;
    public double B01_C_N = MISSING_VALUE;
    public double nN_N = MISSING_VALUE;
    public double F10_C_Cl = MISSING_VALUE;


    public EmbeddedDescriptors(InsilicoMolecule mol) throws GenericFailureException {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws GenericFailureException {

    }

}
