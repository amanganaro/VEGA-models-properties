package insilico.devtox_pg;

import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.devtox_pg.library.PGMoleculeFileSDF;
import insilico.devtox_pg.library.VirtualCompoundLibrary;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import utils.ModelsDeployment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;


public class mainScriptDevToxPG {
    private static final Logger log = LogManager.getLogger(mainScriptDevToxPG.class);

    public static void main(String[] args) throws Exception {

//        VirtualCompoundLibrary.ProcessDirectorySMI(Paths.get("G:\\Drive condivisi\\Tech\\Chemoinformatics\\IRFMN\\In progress\\VEGA - DevTox PG\\PG nuovo\\smiles elaborate e corrette"));
//        InsilicoModel mod = new ismDevToxPG();
//        InsilicoMolecule mol = SmilesMolecule.Convert("CCC=O");
//        System.out.println(mod.Execute(mol).getAssessment());
//        if (1==1) return;

        InsilicoModel model = new ismDevToxPG();
        ModelsDeployment.BuildDataset(model, "out_ts");
        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
        File destinationFile = new File("VegaModels-DevtoxPG\\src\\main\\resources\\data\\ts_devtox_pg.dat");
        try {
            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (Exception ex) {
            log.warn(ex.getMessage());
        }

        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");

        for(String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            for(int i = 0; i < model.GetResultsName().length; i++)
                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
        }


    }

}
