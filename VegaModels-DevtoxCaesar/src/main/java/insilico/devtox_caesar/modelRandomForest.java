package insilico.devtox_caesar;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import weka.classifiers.trees.RandomForest;
import weka.core.Instances;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;


/**
 * Class for the random forest model for Developmental Toxicity model
 */
public class modelRandomForest {

    // data declaration
    private Instances trainingSet;   // training set
    private RandomForest rfTree;     // albero random forest

    public FileWriter fout = null; // crea il test set in formato arff
    public InputStreamReader input3;
    public String buffer = null; // per scrivere il float
    public InputStreamReader input2;


    public modelRandomForest(String arffSource) throws InitFailureException {

        try {

            // legge file del training set
            URL u = null;
            u = getClass().getResource(arffSource);
            InputStreamReader input = new InputStreamReader(u.openStream());

            // costruisce training set
            trainingSet = new Instances(new BufferedReader(input));
            trainingSet.setClassIndex(trainingSet.numAttributes()-1);

            String[] options = new String[8]; // creo un array di stringhe per impostare
                                              // le opzioni delle rf
            options[0] = "-I";
            options[1] = "10";
            options[2] = "-K";
            options[3] = "0";
            options[4] = "-S";
            options[5] = "1";
            options[6] = "-depth";
            options[7] = "12";


            rfTree = new RandomForest();                // costruisco un albero
            rfTree.setOptions(options);                 // setto le opzioni
            rfTree.buildClassifier(trainingSet);        // costruisco il modello

        } catch (Exception e) {
            throw new InitFailureException("Unable to init WEKA model");
        }
    }


    public boolean ExecuteModel(byte[] ARFF_File) throws GenericFailureException {

        ByteArrayInputStream bais = new ByteArrayInputStream(ARFF_File);
        InputStreamReader isr = new InputStreamReader(bais);
        BufferedReader br = new BufferedReader(isr);

        Instances testSet = null;

        double classeAss = -999;
        boolean ReturnVal = false;

        try {
            testSet = new Instances(br);
            testSet.setClassIndex(testSet.numAttributes()-1);
            classeAss = rfTree.classifyInstance(testSet.instance(0)); // il test set contiene sempre una sola istanza

        } catch (Exception e) {
            throw new GenericFailureException("Unable to execute model");
        }

        if(classeAss == 0.0)
            ReturnVal = false;
        else if (classeAss == 1.0)
            ReturnVal = true;

        return ReturnVal;

    }

}
