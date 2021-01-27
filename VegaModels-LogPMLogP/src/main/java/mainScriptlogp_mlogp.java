import insilico.core.model.InsilicoModel;
import insilico.logp_mlogp.ismLogPMLogP;
import lombok.extern.slf4j.Slf4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScriptlogp_mlogp {

    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismLogPMLogP();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-LogPMLogP\\src\\main\\resources\\data\\ts_logp_mlogp.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }

}
