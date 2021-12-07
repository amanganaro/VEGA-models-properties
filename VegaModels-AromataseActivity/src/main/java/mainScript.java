import insilico.aromatase_activity.ismAromataseActivity;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;

public class mainScript {

    public static void main(String[] args) throws Exception {
        ismAromataseActivity mod = new ismAromataseActivity();
        mod.setSkipADandTSLoading(true);
        InsilicoMolecule m = SmilesMolecule.Convert("CCCC(OCC)C(C)");
        InsilicoModelOutput o = mod.Execute(m);
        System.out.println(o.getAssessment());
        for (String s : o.getResults())
            System.out.println(s);
    }
}
