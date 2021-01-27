import insilico.algae_ec50.descriptors.EmbeddedDescriptors;
//import insilico.algae_ec50.ismAlgaeEC50;
import insilico.algae_ec50.ismAlgaeEC50;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import lombok.extern.slf4j.Slf4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScriptAlgaeEC50 {

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismAlgaeEC50();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-AlgaeEC50\\src\\main\\resources\\data\\ts_algae_ec50.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }

//
    }


}
