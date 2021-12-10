package insilico.loael_coral_liver;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import utils.ModelsDeployment;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

public class mainScriptLoaelCoralLiver {

    public static void main(String[] args) throws InitFailureException, GenericFailureException, MalformedURLException, FileNotFoundException {
        InsilicoModel model = new ismLoaelCoralLiver();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-LoaelLiver\\src\\main\\resources\\data\\ts_loael_liver_coral.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
//        ModelsDeployment.TestModelWithTrainingSet(model, "coral_noael_liver_results");

//        ModelsDeployment.TestModelWithTrainingSet(model, "coral_ppb_results");


        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
        smilesList.add("O=C(OC2CC1N(C)C(CC1)C2(C(=O)OC))c3ccccc3");

        for(String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            for(int i = 0; i < model.GetResultsName().length; i++)
                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
        }
    }
}
