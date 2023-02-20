package insilico.vapour_pressure;

import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import utils.ModelsDeployment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

public class mainScriptVapourPressure {

    private static final Logger log = LogManager.getLogger(mainScriptVapourPressure.class);

    public static void main(String[] args) throws Exception {
        InsilicoModel model = new ismVapourPressure();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-VapourPressure\\src\\main\\resources\\data\\ts_vapour_pressure.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }


        List<String> smilesList = new ArrayList<>();
        FileReader fileReader = new FileReader("VegaModels-VapourPressure\\src\\main\\resources\\smiles.txt");
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line = bufferedReader.readLine();
        while (line!=null) {
//            System.out.println(line);
            smilesList.add(line);
            line = bufferedReader.readLine();
        }
//        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
//        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
//        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
//        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");

        for(String smiles : smilesList) {
            System.out.print(smiles);
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++){
//                System.out.println(smiles);
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//            }
        }
    }

}
