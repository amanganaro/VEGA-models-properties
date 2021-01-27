import insilico.daphnia_demetra.ismDaphniaDemetra;
import insilico.core.model.InsilicoModel;
import lombok.extern.slf4j.Slf4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScript {

    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismDaphniaDemetra();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-DaphniaDemetra\\src\\main\\resources\\data\\ts_daphnia_demetra.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }

}
