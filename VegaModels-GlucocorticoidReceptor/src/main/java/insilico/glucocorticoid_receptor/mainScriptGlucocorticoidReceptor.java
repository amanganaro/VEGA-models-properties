package insilico.glucocorticoid_receptor;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.glucocorticoid_receptor.alert.SAGlucocorticoidReceptor;
import insilico.glucocorticoid_receptor.utils.ModelsDeployment;
import lombok.extern.slf4j.Slf4j;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScriptGlucocorticoidReceptor {

    public static void main (String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {
        InsilicoModel model = new ismGlucocorticoidReceptor();
        ModelsDeployment.TestModelWithTrainingSet(model, model.getInfo().getName() + " - New");
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-GlucocorticoidReceptor\\src\\main\\resources\\data\\ts_glucocorticoid_receptor.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex ) {
//            log.warn(ex.getMessage());
//        }

//        SAGlucocorticoidReceptor saGlucocorticoidReceptor = new SAGlucocorticoidReceptor();
//        SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = false;
//        saGlucocorticoidReceptor.Match(SmilesMolecule.Convert("Cc1cc(C)c(c(C)c1)S(=O)(=O)NC(Cn2cc(C#N)c3ccccc23)C(F)(F)F"));
//        return;
    }
}
