import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.report.pdf.ReportPDFSingle;
import insilico.core.model.runner.InsilicoModelWrapper;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.acf.ACFItem;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.daphnia_demetra.ismDaphniaDemetra;
import insilico.core.model.InsilicoModel;
import lombok.extern.slf4j.Slf4j;
import utils.ModelsDeployment;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

@Slf4j
public class mainScript {

    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismDaphniaDemetra();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-DaphniaDemetra\\src\\main\\resources\\data\\ts_daphnia_demetra.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }

        ModelsDeployment.TestModelWithTrainingSet(model, "Daphnia Demetra - New");

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
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.getDescriptorsSize(); i++){
//                System.out.println(model.getDescriptorsNames()[i] + " === " + model.GetDescriptor(i));
//
//            }
//            System.out.println("==============");
//
//        }

//        for(ACFItem acf : model.GetTrainingSet().getACF().getList()){
//            System.out.println(acf.getACF() + " : " + acf.getFrequency());
//        }


    }

}
