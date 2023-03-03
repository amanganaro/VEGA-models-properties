package insilico.carcinogenicity_bb;

import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SABenigniBossa;
import insilico.core.alerts.builders.SABenigniBossaAdditional;
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


public class mainScriptCarcBB {

    private static final Logger log = LogManager.getLogger(mainScriptCarcBB.class);

    
    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismCarcinogenicityBB();
        SABenigniBossa saBenigniBossa = new SABenigniBossa();
        SABenigniBossaAdditional saBenigniBossaAdditional = new SABenigniBossaAdditional();
        for(int i = 0; i < model.GetTrainingSet().getMoleculesSize(); i++ ){
            AlertList result = saBenigniBossa.Calculate(SmilesMolecule.Convert(model.GetTrainingSet().getSMILES(i)));
            AlertList resultAdd = saBenigniBossaAdditional.Calculate(SmilesMolecule.Convert(model.GetTrainingSet().getSMILES(i)));
            System.out.println(model.GetTrainingSet().getSMILES(i));
            if(result.size() > 0)
                System.out.println();
            if(resultAdd.size() > 0)
                System.out.println();
        }
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-CarcinogenicityBB\\src\\main\\resources\\data\\ts_carc_bb.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//
//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
//        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
//        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
//        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
//
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//        }
    }


}
