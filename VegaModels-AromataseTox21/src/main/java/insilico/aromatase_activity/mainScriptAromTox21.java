package insilico.aromatase_activity;

import insilico.aromatase_activity.utils.ModelsDeployment;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;


public class mainScriptAromTox21 {

    private static final Logger log = LogManager.getLogger(mainScriptAromTox21.class);

    
    public static void main (String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {
        InsilicoModel model = new ismAromataseTox21();
//        ModelsDeployment.TestModelWithTrainingSet(model, model.getInfo().getName() + " - New");
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-AromataseTox21\\src\\main\\resources\\data\\ts_aromatase_tox21.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex ) {
            log.warn(ex.getMessage());
        }
    }
}
