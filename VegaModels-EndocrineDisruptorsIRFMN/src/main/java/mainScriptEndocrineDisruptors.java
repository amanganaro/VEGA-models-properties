import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.endocrine_disruptors_irfmn.ismEndocrineDisruptorsIRFMN;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

public class mainScriptEndocrineDisruptors {
    public static void main(String[] args) throws Exception {
        ismEndocrineDisruptorsIRFMN model = new ismEndocrineDisruptorsIRFMN();
//        model.setSkipADandTSLoading(fa  );

        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-EndocrineDisruptorsIRFMN\\src\\main\\resources\\data\\ts_endocrine_disruptors.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }

    }
}
