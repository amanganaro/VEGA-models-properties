package insilico.pparg_up;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.descriptor.blocks.Matrices2D;
import insilico.descriptor.blocks.P_VSA;
import insilico.descriptor.blocks.logP.MLogP;
import insilico.pparg_up.descriptors.EmbeddedDescriptors;
import insilico.pparg_up.utils.ModelsDeployment;
import javassist.runtime.Desc;
import lombok.extern.slf4j.Slf4j;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.nio.file.Files;

@Slf4j
public class mainScript {

    public static void main(String[] args) throws InitFailureException, GenericFailureException, MalformedURLException, FileNotFoundException, DescriptorNotFoundException {

//        double mean = 8.1114599824;
//        double std = 8.4602293961;
//
//        EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert("O=C(Nc1ccc(O)cc1)C"), true);
//        System.out.println("EMBEDDED: " + embeddedDescriptors.MLOGP2);
//        DescriptorBlock block = new MLogP();
//        block.Calculate((SmilesMolecule.Convert("O=C(Nc1ccc(O)cc1)C")));
//        double P_VSA_LogP_4 = Math.pow(block.GetByName("MLogP").getValue(), 2);
//        double scaled_P_VSA_LogP_4 = (P_VSA_LogP_4 - mean)/std;
//        System.out.println("SCALED CALCULATED: " + scaled);

        InsilicoModel model = new ismPPARGup();
//        model.Execute(SmilesMolecule.Convert("O=C(NC1CCCCC1)NS(=O)(=O)c2ccc(cc2)C(=O)C"));


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
//        ModelsDeployment modelsDeployment = new ModelsDeployment().PrintDescriptor(model, "ppargup:descriptors");
        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "ppargup_results_vegaPvsaLogP");
//        modelsDeployment = new ModelsDeployment().PrintScaledDescriptor("ppargup:scaled_descriptors");



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
