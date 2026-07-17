import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.henryslawcoral.ismHenrysLawCoral;

import java.util.ArrayList;
import java.util.List;

public class mainScript {

    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismHenrysLawCoral();

        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
        smilesList.add("O=C(OC2CC1N(C)C(CC1)C2(C(=O)OC))c3ccccc3");

        for(String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            for(int i = 0; i < model.GetResultsName().length; i++) {
                System.out.println(smiles);
                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
            }
        }
    }
}
