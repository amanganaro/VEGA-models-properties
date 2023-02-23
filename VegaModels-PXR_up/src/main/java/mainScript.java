import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.pxr_up.ismPxrUp;
import insilico.pxr_up.utils.ModelsDeployment;

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

        InsilicoModel model = new ismPxrUp();
        InsilicoMolecule mo = SmilesMolecule.Convert("CCC");
        InsilicoModelOutput o= model.Execute(mo);
        if(1==1) return;


//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/ts_pxrup.dat");
//        File destinationFile = new File("VegaModels-PXR_up\\src\\main\\resources\\data\\ts_pxrup.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }
//        new ModelsDeployment().PrintDescriptor(model);


//        new ModelsDeployment().PrintDescriptorBlock(model, new P_VSA());
//        ModelsDeployment modelsDeployment = new ModelsDeployment().PrintDescriptor(model, "pxrup_descriptors");
        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "pxrup_results");

//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("O=C(N(C)C)C"));
//        System.out.println();

//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("CC1=CC=CC=C1OCC(O)C[NH2+]CCOC1=CC=C(C=C1)C(N)=O"));
//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("ClC1=CC=C2N(C3CCN(CCCN4C(=O)NC5=CC=CC=C45)CC3)C(=O)NC2=C1"));
//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("OCC[NH2+]CCNC1=C2C(O)=C3C(O)=CC=C(O)C3=C(O)C2=C(NCC=NCCO)C=C1"));



//        InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert("COC(=O)C1C(C2=CC=CC=C2[N+]([O-])=O)C(C(=O)OC)=C(C)N=C1C"));

//        EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert("CCCCC1=C(C(=O)C2=CC(I)=C(OCC[NH+](CC)CC)C(I)=C2)C2=CC=CC=C2O1"), false);
//        double[] val = embeddedDescriptors.getDescriptors();


//        System.out.println("prediction: " + out.getMainResultValue());
//        System.out.println("prediction^2: " + Math.pow(out.getMainResultValue(), 2));

//        return;
    }
}
