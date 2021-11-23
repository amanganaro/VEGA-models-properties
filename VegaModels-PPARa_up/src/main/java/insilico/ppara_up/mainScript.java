package insilico.ppara_up;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.ppara_up.utils.ModelsDeployment;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScript {

    public static void main(String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {
        InsilicoModel model = new ismPPARAUp();
//        ModelsDeployment.BuildDataset(model, "out_ts");

//        ModelsDeployment modelsDeployment = new ModelsDeployment().PrintDescriptor(model, "descriptors");
        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "results_descriptors_embedded");
    }

}
