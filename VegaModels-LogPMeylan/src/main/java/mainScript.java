import insilico.meylanlogp.ismLogPMeylan;
import insilico.core.model.InsilicoModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import insilico.core.devops.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;


public class mainScript {

    private static final Logger log = LogManager.getLogger(mainScript.class);


    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismLogPMeylan();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-LogPMeylan\\src\\main\\resources\\data\\ts_logp_meylan.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }

}
