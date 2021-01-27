import insilico.core.model.InsilicoModel;
import insilico.daphnia_demetra.ismDaphniaDemetra;
import utils.ModelsDeployment;

import java.util.ArrayList;

public class main {

    public static void main(String[] args) throws Exception {

        ArrayList<InsilicoModel> models = new ArrayList<>();
//
        String[] smilesList = {
                "O=C(OC(C=C)c1ccc2OCOc2(c1))C",
                "O=C2c1ccccc1N(c3c4C=CC(Oc4(cc(OC)c23))(C)C)C",
                "O=[N+]([O-])c1oc(cc1)c2nnc(o2)N"
        };

        models.add(new ismDaphniaDemetra());

        ModelsDeployment.BuildDataset(models, "out_ts");

//        for(InsilicoModel model : models){
//            for(String smiles : smilesList){
//                InsilicoMolecule Mol = SmilesMolecule.Convert(smiles);
//                InsilicoModelOutput out = model.Execute(Mol);
//                System.out.println("RESULTS FOR [" + smiles + "]");
//                for(int i = 0; i < model.GetResultsSize() ; i++){
//                    System.out.println(model.GetResultsName()[i] + " : " + out.getResults()[i]);
//                }
//                System.out.println(out.getMainResultValue());
//            }
//
//        }

    }
}
