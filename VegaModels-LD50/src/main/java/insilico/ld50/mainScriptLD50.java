package insilico.ld50;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.ld50.utils.ModelsDeployment;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScriptLD50 {

    public static void main (String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {

        InsilicoModel model = new ismLD50();
//        ModelsDeployment.BuildDataset(model, "out_ts");
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results", false);
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results_kode_dataset", true);

    }
}
