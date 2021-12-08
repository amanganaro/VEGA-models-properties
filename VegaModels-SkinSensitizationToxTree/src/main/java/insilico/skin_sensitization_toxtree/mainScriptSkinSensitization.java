package insilico.skin_sensitization_toxtree;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.skin_sensitization_toxtree.utils.ModelsDeployment;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScriptSkinSensitization {

    public static void main(String[] args) throws InitFailureException, GenericFailureException, MalformedURLException, FileNotFoundException {

        InsilicoModel model = new ismSkinSensitizationToxTree();

//        ModelsDeployment.BuildDataset(model, "out_ts");

        ModelsDeployment.TestModelWithTrainingSet(model, "results_skin_sensitisation_kode_dataset", true);

//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
//        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
//        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
//        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
//
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++){
//                System.out.println(smiles);
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//            }
//        }
    }
}
