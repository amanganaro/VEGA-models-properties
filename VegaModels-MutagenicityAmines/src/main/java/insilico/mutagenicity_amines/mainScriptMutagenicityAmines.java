package insilico.mutagenicity_amines;

import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.mutagenicity_amines.rules.AromaticAminesSubclassClassifier;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.exception.Intractable;
import utils.ModelsDeployment;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

@Log4j
public class mainScriptMutagenicityAmines {

    public static void main (String[] args) throws InitFailureException, GenericFailureException, FileNotFoundException {
        InsilicoModel model = new ismMutagenicityAmines();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-MutagenicityAmines\\src\\main\\resources\\data\\ts_muta_amines.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }


//        ModelsDeployment.printResultsFroTrainingSet(model);
//
////
        List<String> smilesList = new ArrayList<>();
        smilesList.add("n1cc(nc2c1c3nc(N)n(c3(cc2C))C)c4ccccc4"); // 0 ma è 1
////        smilesList.add("Nc1ccc(cc1)Sc2ccccc2"); // 0 ma è 1
////        smilesList.add("O=S(=O)(O)c2cc(N)ccc2(c1ccc(N)cc1S(=O)(=O)O)"); // -1
////        smilesList.add("O=[N+]([O-])c1ccc(cc1)c3nc2c(ncnc2n3c4ccc(cc4)[N+](=O)[O-])N"); //-1
////
        for(String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            for(int i = 0; i < model.GetResultsName().length; i++)
                System.out.println(smiles + " | " + out.getResults()[i]);
        }

    }
}
