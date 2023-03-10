package insilico.skin_irritation_coral;

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

public class mainScriptSkinIrritationCoral {

    private static final Logger log = LogManager.getLogger(mainScriptSkinIrritationCoral.class);

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismSkinIrritationCoral();
        model.setSkipADandTSLoading(false);
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-SkinIrritationCoral\\src\\main\\resources\\data\\ts_skin_irritation_coral.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }

//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("C=CCCCCCCCCCC");
//        smilesList.add("C(CBr)Cl");
//        smilesList.add("O=C(O)CCCCCCCCC=C");
//        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
//        smilesList.add("O=C(OC2CC1N(C)C(CC1)C2(C(=O)OC))c3ccccc3");
//
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//        }
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.getDescriptorsSize(); i++){
//                System.out.println(model.getDescriptorsNames()[i] + " === " + model.GetDescriptor(i));
//
//            }
//            System.out.println("==============");
//
//        }


    }


}
