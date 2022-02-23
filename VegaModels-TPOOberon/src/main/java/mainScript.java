import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.tpo_oberon.ismTpoOberon;
import insilico.tpo_oberon.utils.ModelsDeployment;
import lombok.extern.slf4j.Slf4j;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

@Slf4j
public class mainScript {

    public static void main(String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException, InterruptedException {
        ismTpoOberon model = new ismTpoOberon();
        model.SetKnnSkipExperimental(true);

//        model.SkipExperimental = true;
//        ModelsDeployment.PrintDescriptor(model, "dataset_descriptors");

//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-TPOOberon\\src\\main\\resources\\data\\ts_tpo_oberon.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex ) {
//            log.warn(ex.getMessage());
//        }

//        Thread.sleep(5000);
        ModelsDeployment.TestModelWithTrainingSet(model, model.getInfo().getName());
//


//        EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert("Oc1ccc(c(O)c1)CCCCCC"));
    }
}
