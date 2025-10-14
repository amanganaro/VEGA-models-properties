package insilico.logk;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;

import insilico.logk.descriptors.EmbeddedDescriptors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;


public class mainScript {
    private static final Logger log = LogManager.getLogger(mainScript.class);

    public static void main(String[] args) throws InitFailureException, GenericFailureException, MalformedURLException, FileNotFoundException {

//        InsilicoModel model = new ismLogK();
////        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("ClC1=C(CCN2CCN(CC2)C2=NSC3=CC=CC=C23)C=C2CC(=O)NC2=C1"));
////        System.out.println();
//
//        ModelsDeployment.TestModelWithTrainingSet(model, "out_ts");
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/ts_logk.dat");
//        File destinationFile = new File("VegaModels-LogK\\src\\main\\resources\\data\\ts_logk.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }



//        ModelsDeployment modelsDeployment = new ModelsDeployment().PrintDescriptor(model, "descriptors");
//        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "results");

//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("COC1=C(OC)C=C2C(N)=NC(=[NH+]C2=C1)N1CCN(CC1)C(=O)C1=CC=CO1"));

//        EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert("CCCCC1=C(C(=O)C2=CC(I)=C(OCC[NH+](CC)CC)C(I)=C2)C2=CC=CC=C2O1"), false);
//        double[] val = embeddedDescriptors.getDescriptors();


//        System.out.println("prediction: " + out.getMainResultValue());
//        System.out.println("prediction^2: " + Math.pow(out.getMainResultValue(), 2));

//        return;
    }
}
