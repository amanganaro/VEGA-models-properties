package insilico.steroidogenesis;

import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
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
//        if (1==1) return;


        iTrainingSet ts = model.GetTrainingSet();
        for (int i=0; i<ts.getMoleculesSize(); i++) {
            String SMI = ts.getSMILES(i);
            if (ts.getMoleculeSet(i) == TrainingSet.MOLECULE_TRAINING)
                model.SetKnnSkipExperimental(true);
            else
                model.SetKnnSkipExperimental(false);
            InsilicoModelOutput o = model.Execute(SmilesMolecule.Convert(SMI));
            System.out.println(SMI + "\t" + o.getMainResultValue());
        }
        if (1==1) return;


        List<String> smilesList = new ArrayList<>();
        URL newData = ismSteroidogenesis.class.getResource("/data/test.txt");
        BufferedReader br = new BufferedReader(new InputStreamReader(newData.openStream()));
        String line;
        while ((line = br.readLine()) != null)
            smilesList.add(line);
        br.close();

        FpKNN k = new FpKNN();
        int n=0;
        for (String smiles : smilesList) {
            double pred = k.Calculate(SmilesMolecule.Convert(smiles), false);
            n++;
//            System.out.println(smiles + "\t" + pred);
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            System.out.println(smiles + "\t" + out.getMainResultValue());
//            for (int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
        }


//        iTrainingSet ts = model.GetTrainingSet();
//        URL KnnData = ismSteroidogenesis.class.getResource("/data/knn_data.txt");
//        KNN knn = new KNN(KnnData);
//
//        for (int i=0; i<ts.getMoleculesSize(); i++) {
//            String SMI = ts.getSMILES(i);
//
//            double[] Descriptors = new double[11];
//
//            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert(SMI));
//            Descriptors[0] = embeddedDescriptors.SRW9;
//            Descriptors[1] = embeddedDescriptors.PubchemFP20;
//            Descriptors[2] = embeddedDescriptors.PubchemFP37;
//            Descriptors[3] = embeddedDescriptors.PubchemFP183;
//            Descriptors[4] = embeddedDescriptors.PubchemFP189;
//            Descriptors[5] = embeddedDescriptors.PubchemFP341;
//            Descriptors[6] = embeddedDescriptors.PubchemFP342;
//            Descriptors[7] = embeddedDescriptors.PubchemFP379;
//            Descriptors[8] = embeddedDescriptors.PubchemFP418;
//            Descriptors[9] = embeddedDescriptors.PubchemFP755;
//
////            if (i==39)
////                System.out.println();
//
//            int prediction = knn.getPredictionCustom(Descriptors, i);
////            int prediction = knn.getPredictionCustom(Descriptors, null);
//            System.out.println(SMI + "\t" + prediction);
//
////            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(SMI));
////            System.out.println(SMI + "\t" + out.getMainResultValue());
//        }
    }
}
