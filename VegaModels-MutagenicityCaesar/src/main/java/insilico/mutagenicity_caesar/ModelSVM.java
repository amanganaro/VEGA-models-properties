package insilico.mutagenicity_caesar;

import insilico.core.exception.InitFailureException;
import insilico.mutagenicity_caesar.libsvm.svm;
import insilico.mutagenicity_caesar.libsvm.svm_model;
import insilico.mutagenicity_caesar.libsvm.svm_node;
import insilico.mutagenicity_caesar.libsvm.svm_parameter;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Class for the SVM model for Mutagenicity
 * 
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class ModelSVM {

    private svm_model model;
    private int predict_probability;
    private double[] ranges;
    
    /**
     * Constructor of the class
     * 
     * @throws InitFailureException
     */
    public ModelSVM(InputStream ModelSource, InputStream ModelRanges) throws InitFailureException {

        try {

            // Inits the model
            model = svm.svm_load_model(ModelSource);
            predict_probability = 0;
            
            // Inits ranges
            String strLine;
            double min, max;
            
            ranges = new double[25];
            BufferedReader br = new BufferedReader(new InputStreamReader(ModelRanges));
            for (int i=0; i<25; i++) {
                strLine = br.readLine();
                String[] buf = strLine.split("\t");
                min = Math.abs(Double.parseDouble(buf[1]));
                max = Math.abs(Double.parseDouble(buf[2]));
                if (min>max)
                    ranges[i] = min;
                else
                    ranges[i] = max;
            }
            
        } catch (Exception e) {
            throw new InitFailureException();
        }
        
    }

    
    public boolean Calculate_Prediction(double[] DescValues) {

        double[] prob_estimates=null;
        double pred;

        int svm_type=svm.svm_get_svm_type(model);
        int nr_class=svm.svm_get_nr_class(model);
        
        svm_node[] x = new svm_node[DescValues.length];
        for (int i=0; i<DescValues.length; i++) {
            x[i] = new svm_node();
            x[i].index = (i+1);
            // Descriptor value is normalized
            x[i].value = (DescValues[i] / ranges[i]);
        }
        
        if (predict_probability==1 && (svm_type== svm_parameter.C_SVC || svm_type==svm_parameter.NU_SVC))
        {
            pred = svm.svm_predict_probability(model,x,prob_estimates);
        } else {
            pred = svm.svm_predict(model,x);
        }

        if (pred==1)
            return false;
        else
            return true;
    }
    
}
