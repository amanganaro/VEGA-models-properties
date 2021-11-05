package insilico.pxr_up.utils;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.pxr_up.descriptors.EmbeddedDescriptors;
import lombok.extern.slf4j.Slf4j;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

@Slf4j
public class ModelsDeployment {

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
            try {
                block.Calculate(SmilesMolecule.Convert(smiles));

            } catch (Exception ex){
                log.warn(ex.getMessage());
            }
        }


    }


    public ModelsDeployment PrintDescriptor(InsilicoModel model, String filename) throws FileNotFoundException {
        List<String> smilesList = new ArrayList<>();
        URL url = (getClass().getResource("/data/PXR_up.csv"));
        StringBuilder stringBuilder = new StringBuilder("#" + "\t" + "Smiles (ionized)");
        String line;
        boolean printHeader = true;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
            while ((line = br.readLine()) != null) {
                if(printHeader){
                    printHeader = false;
                    String[] lineArray = line.split("\t");
                    for(int i = 15; i < lineArray.length; i++) {
                        stringBuilder.append("\t").append(lineArray[i]);
                    }
                } else {
                    smilesList.add(line.split("\t")[2]);
                }
            }

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }

        PrintWriter printWriter = new PrintWriter(filename + ".csv");
        printWriter.print(stringBuilder + "\n");
        printWriter.flush();


        try {
            int index = 1;
            for(String smiles : smilesList) {
                System.out.println("Calculating model descriptors for molecule #" + index + ": " + smiles);


                stringBuilder = new StringBuilder(index + "\t" + smiles);

                EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(SmilesMolecule.Convert(smiles), false);
                for(Double descriptor : embeddedDescriptors.getDescriptors())
                    stringBuilder.append("\t").append(descriptor);
                printWriter.println(stringBuilder);
                printWriter.flush();
                index++;
            }
        } catch (MalformedURLException ex){
            log.warn(ex.getMessage());
        }

        printWriter.flush();
        printWriter.close();




        return this;
    }


    public ModelsDeployment TestModelWithTrainingSet(InsilicoModel model, String filename) throws MalformedURLException, FileNotFoundException, GenericFailureException {
        List<String> smilesList = new ArrayList<>();
        List<String> knimePrediction = new ArrayList<>();
        List<String> knimePrediction_0 = new ArrayList<>();
        List<String> knimePrediction_1 = new ArrayList<>();
        List<String> setType = new ArrayList<>();
        List<String> experimentalValues = new ArrayList<>();
        URL url = (getClass().getResource("/data/PXR_up.csv"));
        String line;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
            br.readLine();
            while ((line = br.readLine()) != null){
                smilesList.add(line.split("\t")[2]);
                setType.add(line.split("\t")[7]);
                experimentalValues.add(line.split("\t")[8]);
                if(line.split("\t")[7].equals("train")){
                    knimePrediction_0.add(line.split("\t")[9]);
                    knimePrediction_1.add(line.split("\t")[10]);
                    if((Double.parseDouble((line.split("\t")[9]))) > (Double.parseDouble((line.split("\t")[10]))))
                        knimePrediction.add("0");
                    else knimePrediction.add("1");

                } else {
                    knimePrediction_0.add(line.split("\t")[13]);
                    knimePrediction_1.add(line.split("\t")[14]);
                    if((Double.parseDouble((line.split("\t")[13]))) > (Double.parseDouble((line.split("\t")[14]))))
                        knimePrediction.add("0");
                    else knimePrediction.add("1");
                }
            }

        } catch (Exception ex){ }
        PrintWriter printWriter = new PrintWriter(filename + ".csv");
        StringBuilder stringBuilder = new StringBuilder("Id\t" + "Smiles\t" + "Experimental Values\t" + "Knime Prediction\t" + "Knime Prediction_0\t" + "Knime Prediction_1\t" + "Vega Prediction\t" + "Vega Prediction_0\t" + "Vega Prediction_1\t" + "Set\n");
//        printWriter.println(stringBuilder);
        int index = 0;
        for(String smiles: smilesList){
            log.info("Calculating PXR_up for molecule #" + (index+1) + ": " + smiles + " ...");
            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
            if(out.getStatus() < 1)
                System.out.println();
            stringBuilder.append(index+1).append("\t")
                    .append(smiles).append("\t")
                    .append(experimentalValues.get(index)).append("\t")
                    .append(knimePrediction.get(index)).append("\t")
                    .append(knimePrediction_0.get(index)).append("\t")
                    .append(knimePrediction_1.get(index)).append("\t")
                    .append(out.getResults()[0]).append("\t")
                    .append(out.getResults()[1]).append("\t")
                    .append(out.getResults()[2]).append("\t")
                    .append(setType.get(index)).append("\n");

            index++;
        }
        printWriter.print(stringBuilder);
        printWriter.close();
        printWriter.flush();

        return this;
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
