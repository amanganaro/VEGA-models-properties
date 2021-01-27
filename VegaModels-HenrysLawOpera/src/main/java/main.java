import insilico.henryslaw.ismHenrysLawOpera;
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
                "COC(=O)C1(O)OCC23C4C(OCC4(C(CC2OC(=O)C(\\\\C)=C\\\\C)OC(C)=O)C(=O)OC)C(O)C(C)(C13)C12OC1(C)C1CC2OC2OC=CC12O",

        };
//
        models.add(new ismHenrysLawOpera());


//        AlgaeDescriptors algaeDescriptors = new AlgaeDescriptors();
//        for(String smiles : smilesList){
//            InsilicoMolecule Mol = SmilesMolecule.Convert(smiles);
//            algaeDescriptors.Calculate(Mol);
//            for(Descriptor d : algaeDescriptors.GetDescriptors()){
//                System.out.println(d.getName() + ": " + d.getValue());
//            }
//        }


//
//        ModelsDeployment.BuildDataset(models, "out_ts");

        for(InsilicoModel model : models){
            for(String smiles : smilesList){
                InsilicoMolecule Mol = SmilesMolecule.Convert(smiles);
                InsilicoModelOutput out = model.Execute(Mol);
                System.out.println("RESULTS FOR [" + smiles + "]");
                if(out.getStatus() == InsilicoModelOutput.OUTPUT_ERROR){
                    System.out.println("FAIL");
                } else {
                    for(int i = 0; i < model.GetResultsSize() ; i++){
                        System.out.println(model.GetResultsName()[i] + " : " + out.getResults()[i]);
                    }
                }
//                System.out.println(out.getMainResultValue());
            }

        }

    }
}
