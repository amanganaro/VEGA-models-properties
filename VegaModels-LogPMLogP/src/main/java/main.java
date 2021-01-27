import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.logp_mlogp.ismLogPMLogP;

public class main {

    public static void main(String[] args) throws Exception {
        String smiles = "O=C1C=CC3(C(=C1)CCC2C4CC(C)C(O)(C(=O)CO)C4(C)(CC(O)C23(F)))(C)";
        InsilicoModel model = new ismLogPMLogP();
        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
        System.out.println(out.getMainResultValue());
    }


}
