import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.endocrine_disruptors_irfmn.ismEndocrineDisruptorsIRFMN;

public class mainScript {
    public static void main(String[] args) throws Exception {
        ismEndocrineDisruptorsIRFMN ed = new ismEndocrineDisruptorsIRFMN();
        ed.setSkipADandTSLoading(true);

        InsilicoMolecule m = SmilesMolecule.Convert("C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)");
        InsilicoModelOutput o = ed.Execute(m);
        System.out.println(o.getResults()[0]);
        System.out.println(o.getResults()[1]);
    }
}
