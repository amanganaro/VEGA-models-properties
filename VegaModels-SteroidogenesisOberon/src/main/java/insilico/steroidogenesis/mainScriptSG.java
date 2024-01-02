package insilico.steroidogenesis;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.knn.insilicoKnnQualitative;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.steroidogenesis.descriptors.EmbeddedDescriptors;
import insilico.steroidogenesis.descriptors.PubchemFingerprinterVega;
import insilico.steroidogenesis.ismSteroidogenesis;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.debug.DebugChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import utils.ModelsDeployment;

import java.io.File;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class mainScriptSG {

    private static final Logger log = LogManager.getLogger(mainScriptSG.class);


    public static void main(String[] args) throws Exception {

        InsilicoModel model = new ismSteroidogenesis();

//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-Steroidogenesis\\src\\main\\resources\\data\\ts_steroidogenesis.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }

//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("CC(C)(C)C1=CC=CC=C1O");
//
//        for (String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            System.out.println("\n"+smiles);
//            for (int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//        }

        iTrainingSet ts = model.GetTrainingSet();
        URL KnnData = ismSteroidogenesis.class.getResource("/data/knn_data_only_train.txt");
        KNN knn = new KNN(KnnData);

        for (int i=0; i<ts.getMoleculesSize(); i++) {
            String SMI = ts.getSMILES(i);

            double[] Descriptors = new double[11];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert(SMI));
            Descriptors[0] = embeddedDescriptors.SRW9;
            Descriptors[1] = embeddedDescriptors.PubchemFP20;
            Descriptors[2] = embeddedDescriptors.PubchemFP37;
            Descriptors[3] = embeddedDescriptors.PubchemFP183;
            Descriptors[4] = embeddedDescriptors.PubchemFP189;
            Descriptors[5] = embeddedDescriptors.PubchemFP341;
            Descriptors[6] = embeddedDescriptors.PubchemFP342;
            Descriptors[7] = embeddedDescriptors.PubchemFP379;
            Descriptors[8] = embeddedDescriptors.PubchemFP418;
            Descriptors[9] = embeddedDescriptors.PubchemFP755;
            System.out.print(SMI + "\t");
            int prediction = knn.getPredictionCustom(Descriptors, i);
//            System.out.println(SMI + "\t" + prediction);

//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(SMI));
//            System.out.println(SMI + "\t" + out.getMainResultValue());
        }
    }
}
