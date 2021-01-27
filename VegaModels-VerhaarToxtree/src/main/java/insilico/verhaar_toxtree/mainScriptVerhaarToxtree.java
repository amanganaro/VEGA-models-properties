package insilico.verhaar_toxtree;

import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import lombok.extern.slf4j.Slf4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScriptVerhaarToxtree {
    public static void main(String[] args) throws Exception {
        InsilicoModel model = new ismVerhaarToxtree();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-VerhaarToxtree\\src\\main\\resources\\data\\ts_verhaar_toxtree.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("CCCC"));
//        ModelsDeployment.FastTestForModel(model);
    }
}
