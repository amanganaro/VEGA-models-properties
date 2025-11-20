package insilico.dio1;

import insilico.core.devops.DeployConfiguration;
import insilico.core.devops.ModelsDeployment;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

public class mainScriptDio1 {

    private static final Logger log = LogManager.getLogger(mainScriptDio1.class);


    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismDio1();

//        ModelsDeployment.BuildDataset(model, "out_ts");
//        System.out.println("done ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-Dio1\\src\\main\\resources\\data\\ts_dio1.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//        if (1==1) return;

        TrainingSet ts = (TrainingSet) model.GetTrainingSet();
        System.out.println("id\tsmiles\tset\texperimental\tpredicted");
        for (int i=0; i<ts.MoleculesSize; i++)
            System.out.println(ts.getId(i) + "\t" + ts.getSMILES(i) + "\t" + ts.getMoleculeSet(i) + "\t"
                    + ts.getExperimentalValue(i) + "\t" + ts.getPredictedValue(i));
        if (1==1) return;


        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=C(O)COc1ccc(cc1)Cl");
        smilesList.add("O=S(=O)(O)c4ccccc4(C=Cc1ccc(cc1)c2ccc(cc2)C=Cc3ccccc3S(=O)(=O)O)");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");

        for(String smiles : smilesList) {
            System.out.println("* " + smiles);
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            for(int i = 0; i < model.GetResultsName().length; i++){
                System.out.println("   " + model.GetResultsName()[i] + " | " + out.getResults()[i]);
            }
        }


    }
}
