package insilico.ahr_up;

import insilico.ahr_up.utils.ModelsDeployment;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScript {

    public static void main(String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {

        InsilicoModel model = new ismAhrUp();
//        ModelsDeployment.BuildDataset(model, "out_ts");

        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "results_ahrup_embedded_descriptors");

    }
}
