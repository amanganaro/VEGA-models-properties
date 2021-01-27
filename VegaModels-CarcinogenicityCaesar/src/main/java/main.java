import insilico.carcinogenicity_caesar.descriptors.EmbeddedDescriptors;
import insilico.core.model.InsilicoModel;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;

import java.util.ArrayList;

public class main {

    public static void main(String[] args) throws Exception {

        ArrayList<InsilicoModel> models = new ArrayList<>();
////
        String[] smilesList = {
                "O=[N+]([O-])c1oc(cc1)c2cnc(NN=C(C)C)s2",
                "N#CC",
                "ON=C(C)C",
                "O=C(OC(C=C)c1ccc2OCOc2(c1))C"

        };
//
//        models.add(new ismCarcinogenicityCaesar());

        for(String smiles : smilesList) {
            InsilicoMolecule Mol = SmilesMolecule.Convert(smiles);
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();
            embeddedDescriptors.CalculateAllDescriptors(Mol);
            System.out.println("[" + smiles + "]");
            System.out.println("MATS2P : " + embeddedDescriptors.MATS2p);
            System.out.println("D/Dr6 : "+ embeddedDescriptors.DDr6);
            System.out.println("EEig10ed : "+ embeddedDescriptors.EEig10ed);
            System.out.println("ESpm11ed : "+ embeddedDescriptors.ESpm11ed);
            System.out.println("ESpm9dm :"+ embeddedDescriptors.ESpm9dm);

            System.out.println("GGI2 : "+ embeddedDescriptors.GGI2);
            System.out.println("JGI6 : "+ embeddedDescriptors.JGI6);
            System.out.println("PW5 : "+ embeddedDescriptors.PW5);
        }



//
//        ModelsDeployment.BuildDataset(models, "out_ts");
//
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
////
//        }


    }

}
