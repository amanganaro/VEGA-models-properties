package insilico.logp_alogp;

import insilico.core.devops.ModelsDeployment;
import insilico.core.model.InsilicoModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;


public class mainScriptALogP {
    private static final Logger log = LogManager.getLogger(mainScriptALogP.class);

    public static void main(String[] args) throws Exception {
        InsilicoModel model = new ismLogPALogP();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-LogPALogP\\src\\main\\resources\\data\\ts_logp_alogp.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
//        ModelsDeployment.FastTestModel(model);
    }
}
