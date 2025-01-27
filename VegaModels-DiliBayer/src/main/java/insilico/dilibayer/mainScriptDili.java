package insilico.dilibayer;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.python.CdddDescriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import utils.ModelsDeployment;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

public class mainScriptDili {

    private static final Logger log = LogManager.getLogger(mainScriptDili.class);


    public static void main(String[] args) throws GenericFailureException, InitFailureException, IOException, URISyntaxException, InterruptedException {
        InsilicoModel model = new ismDiliBayer(true, null);

//        iTrainingSet ist = model.GetTrainingSet();
//        for(int i=0; i < ist.getMoleculesSize(); i++){
//            System.out.println(ist.getSMILES(i)+" "+ist.getPredictedValueFormatted(i));
//        }
//
//        if(1==1)
//            return;
//        ModelsDeployment.printResultsFroTrainingSet(model);
//        if(1==1)
//            return;
//
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-DiliBayer\\src\\main\\resources\\data\\ts_dili_bayer.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//        if(1==1)
//            return;
//
//        model.setSkipADandTSLoading(true);

        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");

        CdddDescriptors cdddDescriptors=null;

        if(InsilicoModelPython.class.isAssignableFrom(model.getClass())){
             cdddDescriptors = new CdddDescriptors(smilesList, false);
            ((ismDiliBayer) model).setDescriptorGenerator(cdddDescriptors);
            boolean descriptorOK = cdddDescriptors.calculateDescriptors();
        }

        for (String smiles : smilesList) {
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            System.out.println("\n"+smiles);
            System.out.println("\n"+out.getAssessment() + "\n");
            System.out.println("\n"+out.getAssessmentVerbose() + "\n");
            for (int i = 0; i < model.GetResultsName().length; i++)
                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
            for (int i = 0; i < model.GetADItemsName().length; i++)
                System.out.println(model.GetADItemsName()[i] + " | " + out.getADIndex().get(i).GetIndexValueFormatted() );
        }

        cdddDescriptors.dispose();
    }
}
