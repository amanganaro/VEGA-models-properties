package insilico.skin_irritation.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.Rings;
import insilico.core.descriptor.blocks.weights.basic.WeightsElectronegativity;
import insilico.core.descriptor.blocks.weights.basic.WeightsIonizationPotential;
import insilico.core.descriptor.blocks.weights.basic.WeightsMass;
import insilico.core.descriptor.blocks.weights.other.WeightsHydrophobicityGC;
import insilico.core.descriptor.blocks.weights.other.WeightsVertexDegree;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.GhoseCrippenACF;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
//import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class EmbeddedDescriptors {

    private static final Logger log = LogManager.getLogger(EmbeddedDescriptors.class);

    private int MISSING_VALUE = -999;

    public double SM1_Dzi = MISSING_VALUE;
    public double Eig09_EAri = MISSING_VALUE;
    public double NdsCH = MISSING_VALUE;
    public double X3Av = MISSING_VALUE;
    public double NssS = MISSING_VALUE;
    public double nRNR2 = MISSING_VALUE;
    public double CATS2D_01_DA = MISSING_VALUE;
    public double IC3 = MISSING_VALUE;
    public double BB_SA44 = MISSING_VALUE;
    public double nArCHO = MISSING_VALUE;
    public double nRCO = MISSING_VALUE;
    public double SM03_AEAdm = MISSING_VALUE;
    public double B03_N_O = MISSING_VALUE;
    public double B04_C_N = MISSING_VALUE;
    public double Eig05_EAdm = MISSING_VALUE;
    public double C_026 = MISSING_VALUE;
    public double nCrs = MISSING_VALUE;
    public double SpMaxA_D_Dt = MISSING_VALUE;
    public double H_052 = MISSING_VALUE;
    public double nN = MISSING_VALUE;
    public double CATS2D_09_AL = MISSING_VALUE;
    public double Eig12_EAdm = MISSING_VALUE;
    public double Eig15_EAdm = MISSING_VALUE;
    public double CATS2D_08_NL = MISSING_VALUE;
    public double SpMin6_Bhs = MISSING_VALUE;
    public double SpMin6_Bhi = MISSING_VALUE;
    public double piPC10 = MISSING_VALUE;
    public double B06_O_O = MISSING_VALUE;
    public double Eig09_AEAri = MISSING_VALUE;
    public double SpMAD_AEAdm = MISSING_VALUE;
    public double B05_O_Cl = MISSING_VALUE;
    public double nR5 = MISSING_VALUE;
    public double nArOH = MISSING_VALUE;
    public double B04_Cl_Cl = MISSING_VALUE;
    public double IC1 = MISSING_VALUE;
    public double MPC09 = MISSING_VALUE;
    public double B09_O_O = MISSING_VALUE;
    public double CATS2D_08_PL = MISSING_VALUE;
    public double F08_C_N = MISSING_VALUE;
    public double GATS5s = MISSING_VALUE;
    public double F03_C_Cl = MISSING_VALUE;
    public double H_051 = MISSING_VALUE;
    public double B06_C_N = MISSING_VALUE;
    public double GATS1e = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule mol) throws GenericFailureException {
        CalculateAllDescriptors(mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule mol) throws GenericFailureException {

    }

}
