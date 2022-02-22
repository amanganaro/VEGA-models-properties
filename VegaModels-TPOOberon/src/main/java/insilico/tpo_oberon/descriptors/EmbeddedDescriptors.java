package insilico.tpo_oberon.descriptors;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.descriptor.blocks.*;
import insilico.descriptor.blocks.logP.MLogP;

public class EmbeddedDescriptors {

    private final double MISSING_VALUE = -999;

    public double GATS1e;
    public double nArOH;
    public double CATS2D_02_DL;
    public double MATS1e;
    public double MATS1s;
    public double C_026;
    public double CATS2D_03_DL;
    public double B10_C_C;
    public double MATS1p;
    public double nCb_;
    public double nX;
    public double Uc;
    public double P_VSA_i_1;
    public double SpMAD_B_v_;
    public double nCbH;
    public double GATS1s;
    public double MATS1m;
    public double MLOGP;
    public double SpMax2_Bh_s;
    public double Eta_C_A;

    public double[] descriptorsArray;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws DescriptorNotFoundException {
        CalculateDescriptors(mol);
    }

    private void CalculateDescriptors(InsilicoMolecule mol) throws DescriptorNotFoundException {

        DescriptorBlock block = new AutoCorrelation();
        block.Calculate(mol);
        GATS1e = block.GetByName("GATS1e").getValue();
        MATS1e = block.GetByName("MATS1e").getValue();
        MATS1m = block.GetByName("MATS1m").getValue();
        MATS1p = block.GetByName("MATS1p").getValue();
        MATS1s = block.GetByName("MATS1s").getValue();
        GATS1s = block.GetByName("GATS1s").getValue();


        block = new FunctionalGroups();
        block.Calculate(mol);
        nArOH = block.GetByName("nArOH").getValue();
        nCb_ = block.GetByName("nCb–").getValue();
        nCbH = block.GetByName("nCbH").getValue();


        block = new Cats2D();
        block.Calculate(mol);
        CATS2D_02_DL = block.GetByName("CATS2D_02_DL").getValue();
        CATS2D_03_DL = block.GetByName("CATS2D_03_DL").getValue();

        block = new AtomCenteredFragments();
        block.Calculate(mol);
        C_026 = block.GetByName("C-026").getValue();

        block = new AtomPairs2D();
        block.Calculate(mol);
        B10_C_C = block.GetByName("B10[C-C]").getValue();

        block = new Constitutional();
        block.Calculate(mol);
        nX = block.GetByName("nX").getValue();

        block = new MolecularProperties();
        block.Calculate(mol);
        Uc = block.GetByName("Uc").getValue();

        block = new P_VSA();
        block.Calculate(mol);
        P_VSA_i_1 = block.GetByName("P_VSA_i_1").getValue();

        block = new Matrices2D();
        block.Calculate(mol);
        SpMAD_B_v_ = block.GetByName("SpMAD_B(v)").getValue();


        block = new BurdenEigenvalues();
        block.Calculate(mol);
        SpMax2_Bh_s = block.GetByName("SpMax2_Bh(s)").getValue();

        block = new MLogP();
        block.Calculate(mol);
        MLOGP = block.GetByName("MLOGP").getValue();

        block = new EtaIndices();
        block.Calculate(mol);
        Eta_C_A = block.GetByName("Eta_C_A").getValue();

//        getETA(mol);
        descriptorsArray = new double[]{
                GATS1e,
                nArOH,
                CATS2D_02_DL,
                MATS1e,
                MATS1s,
                C_026,
                CATS2D_03_DL,
                B10_C_C,
                MATS1p,
                nCb_,
                nX,
                Uc,
                P_VSA_i_1,
                SpMAD_B_v_,
                nCbH,
                GATS1s,
                MATS1m,
                MLOGP,
                SpMax2_Bh_s,
                Eta_C_A
        };

    }

    
}
