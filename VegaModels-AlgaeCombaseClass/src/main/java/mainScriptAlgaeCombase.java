import insilico.algae_combaseclass.ismAlgaeCombaseClass;
import insilico.core.localization.StringSelectorCore;
import insilico.core.model.InsilicoModel;
import lombok.extern.log4j.Log4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Log4j
public class mainScriptAlgaeCombase {

    public static void main(String[] args) throws Exception {

        StringSelectorCore.SetLanguage("it");
        InsilicoModel model = new ismAlgaeCombaseClass();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-AlgaeCombaseClass\\src\\main\\resources\\data\\ts_algae_combaseClass.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }
    }
}
