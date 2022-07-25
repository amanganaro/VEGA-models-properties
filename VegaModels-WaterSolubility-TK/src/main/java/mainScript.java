import insilico.core.model.InsilicoModel;
//import utils.ModelsDeployment;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.watersolubility_tk.ismWaterSolubilityIRFMN;
import lombok.extern.log4j.Log4j;

import java.util.ArrayList;
import java.util.List;

@Log4j
public class mainScript {

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismWaterSolubilityIRFMN();
//        ModelsDeployment.TestModelWithTrainingSet(model, "Water solubility model (IRFMN) - New");
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-WaterSolubility\\src\\main\\resources\\data\\ts_watersolubility_irfmn.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//
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

//        for (int i=0; i<model.GetTrainingSet().getMoleculesSize(); i++) {
//            String SMI = model.GetTrainingSet().getSMILES(i);
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(SMI));
//            System.out.println(SMI+"\t" + model.GetTrainingSet().getPredictedValue(i) + "\t" + out.getMainResultValue());
//        }

    }

}
