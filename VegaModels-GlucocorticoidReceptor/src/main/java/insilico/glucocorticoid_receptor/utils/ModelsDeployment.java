package insilico.glucocorticoid_receptor.utils;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;


public class ModelsDeployment {
    private static final Logger log = LogManager.getLogger(ModelsDeployment.class);

    public void PrintDescriptorBlock(InsilicoModel model, DescriptorBlock block){
        List<String> smilesList = new ArrayList<>();
        URL url = (getClass().getResource("/data/SQfu.csv"));
        StringBuilder stringBuilder = new StringBuilder("#" + "\t" + "Smiles (ionized)");
        String line;
        boolean printHeader = true;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
            while ((line = br.readLine()) != null) {
                if(printHeader){
                    printHeader = false;
                    String[] lineArray = line.split("\t");
                    for(int i = 7; i < lineArray.length; i++) {
                        stringBuilder.append("\t").append(lineArray[i]);
                    }
                } else {
                    smilesList.add(line.split("\t")[2]);
                }
            }

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }


        for(String smiles: smilesList){
            block.Calculate(SmilesMolecule.Convert(smiles));
        }


    }

    public static void PrintDescriptor(InsilicoModel model, String filename) throws FileNotFoundException {
        List<String> smilesList = new ArrayList<>();


        String datasetUrl = model.getInfo().getTrainingSetURL().split("\\.")[0] + ".txt";
        URL url = ModelsDeployment.class.getResource(datasetUrl);

        StringBuilder stringBuilder = new StringBuilder("#" + "\t" + "Smiles (ionized)\t");
        for(String descriptorName : model.getDescriptorsNames())
            stringBuilder.append(descriptorName).append("\t");


        PrintWriter printWriter = new PrintWriter(filename + ".csv");
        printWriter.print(stringBuilder + "\n");
        printWriter.flush();


        String line;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
            br.readLine();
            while ((line = br.readLine()) != null) {
                smilesList.add(line.split("\t")[2]);
            }

            for(String smiles: smilesList){
                InsilicoModelOutput curOutput = model.Execute(SmilesMolecule.Convert(smiles));
                stringBuilder = new StringBuilder(smilesList.indexOf(smiles) + 1).append("\t").append(smiles);
                System.out.println("Printing Descriptors for #:" + (smilesList.indexOf(smiles) + 1) + " - " + model.getInfo().getKey() + " - " + smiles);
                for(int i = 0; i < model.getDescriptorsSize(); i++)
                    stringBuilder.append("\t").append(model.GetDescriptor(i));

                printWriter.print(stringBuilder + "\n");
                printWriter.flush();
            }

            printWriter.flush();
            printWriter.close();

        } catch (GenericFailureException | IOException e) {
            log.warn(e.getClass() + ": " + e.getMessage());
        }


    }


    public static void TestModelWithTrainingSet(InsilicoModel model, String csvFilename) throws MalformedURLException, FileNotFoundException, GenericFailureException {
        List<String> smilesList = new ArrayList<>();

        String datasetUrl = model.getInfo().getTrainingSetURL().split("\\.")[0] + ".txt";
        URL url = ModelsDeployment.class.getResource(datasetUrl);
        String line;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
            br.readLine();
            while ((line = br.readLine()) != null){
                smilesList.add(line.split("\t")[2]);
//                knimePredictionList.add(line.split("\t")[4]);
            }

        } catch (Exception ex){ }
        PrintWriter printWriter = new PrintWriter(csvFilename +  ".csv");
        StringBuilder stringBuilder = new StringBuilder("Smiles"  + "\t");

        for(String resultName : model.GetResultsName())
            stringBuilder.append(resultName).append("\t");

        stringBuilder.append("\n");
        printWriter.print(stringBuilder);
        printWriter.flush();
//        printWriter.println(stringBuilder);
        int index = 0;
        for(String smiles: smilesList){
            log.info((index+1) + ": " + smiles + " ...");
            SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = false;
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            if(index == 348)
//                System.out.println();
            stringBuilder = new StringBuilder(smiles).append("\t");
            if(out.getStatus() > -1){
                for(String results : out.getResults())
                    stringBuilder.append(results).append("\t");
            }
            else stringBuilder.append(out.getErrMessage());



            stringBuilder.append("\n");

            printWriter.print(stringBuilder);
            printWriter.flush();
            index++;
        }


        printWriter.print(stringBuilder);
        printWriter.flush();
        printWriter.close();
    }



    /**
     * Process all available models to generate TS molecule images.
     * List of all models can be retrieved from ModelsList class.
     *
     * @param Models list of models to be processed
     * @param DestDir destination folder to save the files
     * @throws Exception
     */
    public static void BuildPNG(ArrayList<InsilicoModel> Models, String DestDir) {

//        int index = 0;
//        for (InsilicoModel model : Models){
//            System.out.println("Processing model no. " + ++index + ": " + model.getInfo().getDescription() + " (" + model.getInfo().getKey() + ")");
//            try {
//                String path = DestDir + model.getInfo().getTrainingSetPngURL();
//
//                new File(path).mkdirs();
//                Depiction.SaveTSMoleculesAsPNG(model.GetTrainingSet(), path);
//
//            } catch (Throwable e){
//                System.out.println("=== ERROR while processing model - " + e.getMessage());
//            }
//        }
//        System.out.println();
//        System.out.println(index + " model processed");
    }


    /**
     * Process all available models to prepare the trainingsets.
     * List of all models can be retrieved from ModelsList class.
     *
     * @param Models list of models to be processed
     * @param DestDir destination folder to move the files (null to leave all
     * files in current directory)
     * @throws Exception
     */
    public static void BuildDataset(ArrayList<InsilicoModel> Models, String DestDir) throws Exception {
       int index = 0, errCount = 0;
       for(InsilicoModel model : Models) {
           int[] indexAndErrors = ProcessModel(model, DestDir);
           index += indexAndErrors[0];
           errCount += indexAndErrors[1];
       }

        System.out.println();
        System.out.println(index + " model processed");
        System.out.println(errCount + " model with errors");

    }

    public static void BuildDataset(InsilicoModel model, String DestDir) {

        int[] indexAndErrors = ProcessModel(model, DestDir);
        int index = indexAndErrors[0], errCount = indexAndErrors[1];


        System.out.println();
        System.out.println(index + " model processed");
        System.out.println(errCount + " model with errors");
    }


    public static void ListAllDescriptors() throws Exception {
//        ArrayList<InsilicoModel> models = ModelsList.getAllModels();
//        for (InsilicoModel m : models) {
//            System.out.println(m.getInfo().getName());
//            for (String name : m.getDescriptorsNames())
//                System.out.println(name);
//            System.out.println();
//        }
    }

    public static void FastTestModel(InsilicoModel model){

        List<String> smilesList = new ArrayList<>();
        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
        smilesList.add("O=C(OC2CC1N(C)C(CC1)C2(C(=O)OC))c3ccccc3");

        try {

            System.out.println("== " + model.getInfo().getName() + " ==");
            System.out.println();
            for(String smiles : smilesList) {
                System.out.println("== [" + smiles + "] ==");
                InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
                System.out.println("---> Model results");
                for(int i = 0; i < model.GetResultsName().length; i++)
                    System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
                System.out.println("---> Descriptors");
                for(int i = 0; i < model.getDescriptorsSize(); i++){
                    System.out.println(model.getDescriptorsNames()[i] + " === " + model.GetDescriptor(i));

                }
                System.out.println("==============");
            }
        } catch (GenericFailureException exception) {
            System.out.println(exception.getMessage());
        }

    }

    private static int[] ProcessModel(InsilicoModel model, String DestDir) {
        int[] result = new int[2];

        // result[0] index, result[1] errCount

        System.out.println("* Processing model no. " + ++result[0] + ": " + model.getInfo().getName() + " (" + model.getInfo().getKey() + ")...");
        try {
            // TrainingSet generation
            model.ProcessTrainingSet();

            // move files
            if (DestDir != null){
                String[] buf = model.getInfo().getTrainingSetURL().split("/");
                String FileName = buf[buf.length -1];
                StringBuilder FinalDir = new StringBuilder(DestDir);
                for (String s : buf)
                    if (!s.isEmpty())
                        FinalDir.append("/").append(s);
                System.out.println(" === Creating directory: " + FinalDir);
                new File(String.valueOf(FinalDir)).mkdirs();
                Files.move(Paths.get(FileName), Paths.get(FinalDir + "/" + FileName), StandardCopyOption.REPLACE_EXISTING);
            }

        } catch (Throwable e) {
            result[1]++;
            System.out.println("*** ERROR while processing model - " + e.getMessage());
        }

        return result;

    }


}
