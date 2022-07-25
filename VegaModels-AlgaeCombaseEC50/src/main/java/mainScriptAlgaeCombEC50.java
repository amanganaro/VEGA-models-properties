import insilico.algae_combaseEC50.ismAlgaeCombaseEC50;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import lombok.extern.log4j.Log4j;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Log4j
public class mainScriptAlgaeCombEC50 {

    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismAlgaeCombaseEC50();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/"  + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-AlgaeCombaseEC50\\src\\main\\resources\\data\\ts_algae_combaseEC50.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
//        }

//        String smiles = "O=C1C=CC3(C(=C1)CCC2C4CC(C)C(O)(C(=O)CO)C4(C)(CC(O)C23(F)))(C)";
//        InsilicoModel model = new ismAlgaeCombaseEC50();
//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//        System.out.println(out.getMainResultValue());
    }
    }


}
