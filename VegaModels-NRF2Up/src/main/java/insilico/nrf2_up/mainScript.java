package insilico.nrf2_up;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.nrf2_up.utils.ModelsDeployment;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScript {

    public static void main(String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {

        ismNRF2Up model = new ismNRF2Up();
        InsilicoMolecule mo = SmilesMolecule.Convert("CCC");
        InsilicoModelOutput o= model.Execute(mo);
        if(1==1) return;

        model.setLoadDescriptorsFromFile(false);

//        ModelsDeployment.BuildDataset(model, "out_ts");

        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "results_nrf2up_descriptors_embedded");

//        model.setLoadDescriptorsFromFile(false);
//        modelsDeployment.TestModelWithTrainingSet(model, "results_nrf2up_descriptor_block");


    }
}
