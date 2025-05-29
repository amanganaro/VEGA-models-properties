package insilico.ld50quail;

import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;


public class mainScriptLD50quail {
    private static final Logger log = LogManager.getLogger(mainScriptLD50quail.class);

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismLD50Quail();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-LD50_Quail\\src\\main\\resources\\data\\ts_ld50quail.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }

        List<String> smilesList = new ArrayList<>();
//        smilesList.add("N#CC(=NOP(=S)(OCC)OCC)c1ccccc1Cl");
//        smilesList.add("C(=NC)=S");
        for (int i=0; i<model.GetTrainingSet().getMoleculesSize(); i++)
            smilesList.add(model.GetTrainingSet().getSMILES(i));

        for(String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            System.out.println(smiles + "\t" + out.getResults()[0]);
//            for(int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
        }

    }


}
