import insilico.carcinogenicity_antares.ismCarcinogenicityAntares;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;

import java.util.ArrayList;

public class main {

    public static void main(String[] args) throws Exception {
        ArrayList<InsilicoModel> models = new ArrayList<>();
////
        String[] smilesList = {
//                "O(CC(C)C)CC(CN(c1ccccc1)Cc2ccccc2)[NH+]3CCCC3",
//                "O=C([O-])C=C",
//                "Nc1ccc(cc1)C(c2ccc(N)cc2)=C3C=CC(=[NH2+])C=C3",
                "O=[N+]([O-])c1oc(cc1)c2cnc(NN=C(C)C)s2"
        };
//
        models.add(new ismCarcinogenicityAntares());


//
//        ModelsDeployment.BuildDataset(models, "out_ts");
//
        for(InsilicoModel model : models){
            for(String smiles : smilesList){
                InsilicoMolecule Mol = SmilesMolecule.Convert(smiles);
                InsilicoModelOutput out = model.Execute(Mol);
                System.out.println("RESULTS FOR [" + smiles + "]");
                for(int i = 0; i < model.GetResultsSize() ; i++){
                    System.out.println(model.GetResultsName()[i] + " : " + out.getResults()[i]);
                }
                System.out.println(out.getMainResultValue());
            }
//
        }

    }
}

