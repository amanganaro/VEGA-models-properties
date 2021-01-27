package utils;

import insilico.core.model.InsilicoModel;
//import ModelsList;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

public class ModelsDeployment {

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
