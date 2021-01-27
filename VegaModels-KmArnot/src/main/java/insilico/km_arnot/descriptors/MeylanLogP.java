package insilico.km_arnot.descriptors;

import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAMeylanLogPAdditionalFragments;
import insilico.core.alerts.builders.SAMeylanLogPCorrectionFragments;
import insilico.core.alerts.builders.SAMeylanLogPFragments;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.old.AtomCenteredFragments;
import insilico.core.descriptor.blocks.old.weight.GhoseCrippenWeights;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Meylan LogP descriptor (from Kowwin).
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class MeylanLogP extends DescriptorBlock {

    private final static boolean FRAGS_TO_OUT = false;

    private final static long serialVersionUID = 1L;
    private final static String BlockName = "Meylan LogP (Kowwin)";

    /**
     * Constructor. This should not be used, no weight is specified. The
     * overloaded constructors should be used instead.
     */
    public MeylanLogP() {
        super();
        this.Name = MeylanLogP.BlockName;
    }



    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        Add("LogP", "Meylan LogP (Kowwin)");
        SetAllValues(Descriptor.MISSING_VALUE);
    }


    /**
     * Calculate descriptors for the given molecule.
     *
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {
        this.Calculate(mol, null);
    }


    /**
     * Overload of Calculate(), to be used if needed matrix is passed (to
     * fasten calculation)
     *
     * @param mol molecule to be calculated
     * @param ACF Atom Centered Fragments descriptor block for the current
     * molecule
     */
    public void Calculate(InsilicoMolecule mol, DescriptorBlock ACF) {

        // Generate/clears descriptors
        GenerateDescriptors();

        // Calculates fragments
        try {
            IAtomContainer m = mol.GetStructure();

            SAMeylanLogPFragments Frags = new SAMeylanLogPFragments();
            AlertList f_base = Frags.Calculate(mol);
            SAMeylanLogPAdditionalFragments FragsAdd = new SAMeylanLogPAdditionalFragments();
            AlertList f_add = FragsAdd.Calculate(mol);
            SAMeylanLogPCorrectionFragments FragsCorr = new SAMeylanLogPCorrectionFragments();
            AlertList f_corr = FragsCorr.Calculate(mol);

            if (FRAGS_TO_OUT) {
                for (Alert a : f_base.getSAList())
                    System.out.println(a.getId() + "\t" + a.getName() + "\t" + a.getDescription());
                for (Alert a : f_add.getSAList())
                    System.out.println(a.getId() + "\t" + a.getName() + "\t" + a.getDescription());
                for (Alert a : f_corr.getSAList())
                    System.out.println(a.getId() + "\t" + a.getName() + "\t" + a.getDescription());
            }

            double Coefficient = Frags.GetCoefficient() + FragsAdd.GetCoefficient();
            double Correction = FragsCorr.GetCoefficient();

            // main eq
            double Meylanlogp = 0.2290 + Coefficient + Correction;

            // lower bound
            if (Meylanlogp < -5.0)
                Meylanlogp = -5.0;

            SetByName("LogP", Meylanlogp);

        } catch (InitFailureException | InvalidMoleculeException | GenericFailureException e) {
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF;
        if (ACF != null)
            CurACF = ACF;
        else {
            CurACF = new AtomCenteredFragments();
            CurACF.Calculate(mol);
        }

        double LogP = 0;
        double MR = 0;
        double[] Frags = CurACF.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
            MR += Frags[i] * GhoseCrippenWeights.GetMolarRefractivity(i);
        }

        SetByName("ALogP", LogP);
        SetByName("AMR", MR);
    }


    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    @Override
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        MeylanLogP block = new MeylanLogP();
        block.CloneDetailsFrom(this);
        return block;
    }


}
