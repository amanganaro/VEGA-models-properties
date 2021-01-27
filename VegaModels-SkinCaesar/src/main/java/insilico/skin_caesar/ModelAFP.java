package insilico.skin_caesar;

import insilico.core.exception.InitFailureException;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;


/**
 * Class for the AFP model for Skin Sensitization
 * 
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class ModelAFP {

    private final static double Th_Accept = 0.5;        // Acceptance threshold
    private final static double Th_Exclusion = 0.001;   // Exclusion threshold
    
    private AFPRule[] Rules;
    private int Rules_Num;
   
    
    /**
     * Constructor of the class
     * 
     * @throws InitFailureException
     */
    public ModelAFP(String RulesFile) throws InitFailureException {

        try {

            // Loads rules for the model
            LoadRules(RulesFile);

            // Loads training set
            // TO DO
            
        } catch (Exception e) {
            throw new InitFailureException();
        }
        
    }
    
    
    /**
     * Loads the rules for the AFP from the csv file
     * 
     * @throws Exception
     */
    private void LoadRules(String RulesFile) throws Exception {

	int n=0, i=0, j=0, k=0;
        AFPRule[] BufArray;
        int Idx=0;
        
        Rules = new AFPRule[0];

        URL u;

        u = getClass().getResource(RulesFile);
        InputStream is = u.openStream();
        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String strLine;

        // Skips first line (headers)
        strLine = br.readLine();

        while ((strLine = br.readLine()) != null)   {

            // Gets number of tokens
            String[] buf = strLine.split(",");

            // Resizes buffer array for the new line to be added
            BufArray = new AFPRule[Idx+1];
            for (i=0;i<Idx;i++)
                BufArray[i] = Rules[i];

            // Adds the new line read
            BufArray[Idx] = new AFPRule();
            j=1;
            double O_Inact = Double.parseDouble(buf[j++]);
            double O_Act = Double.parseDouble(buf[j++]);
            BufArray[Idx].SetRule(Idx, O_Act, O_Inact);
            for (k=0; k< AFPRule.DescriptorsNum; k++) {
                double A = Double.parseDouble(buf[3 + k * 8]);
                double B = Double.parseDouble(buf[3 + k * 8 + 1]);
                double C = Double.parseDouble(buf[3 + k * 8 + 2]);
                double D = Double.parseDouble(buf[3 + k * 8 + 3]);
                double yA = Double.parseDouble(buf[3 + k * 8 + 4]);
                double yB = Double.parseDouble(buf[3 + k * 8 + 5]);
                double yC = Double.parseDouble(buf[3 + k * 8 + 6]);
                double yD = Double.parseDouble(buf[3 + k * 8 + 7]);
                if ( (A==-999)||(B==-999)||(C==-999)||(D==-999) ) {
                    BufArray[Idx].Used[k] = false;
                    BufArray[Idx].AddDescriptorRule(k, -999, -999, -999, -999,
                        -999, -999, -999, -999);
                } else {
                    BufArray[Idx].Used[k] = true;
                    BufArray[Idx].AddDescriptorRule(k, A, B, C, D, yA, yB, yC, yD);
                }
            }

            Rules = BufArray;
            Idx++;
        }

        Rules_Num = Idx;
    }
    
    
    /**
     * Calculates the prediction for a given molecule, passed as a set of
     * descriptors
     * 
     * @param DescValues    array of double containing descriptors
     * @return  the prediction as a ModelAFPResult object
     */
    public ModelAFPResult Calculate_Prediction(double[] DescValues) {
        
        double[] u_results;
        double[] O_Act;
        double[] O_Inact;
        double min_u = 0, cur_u = 0;
        double Sum_Act = 0,Sum_Inact = 0;
        double Sum_u = 0;
        double Final_O_Act,Final_O_Inact;
        
        ModelAFPResult Result = null;
        
        u_results = new double[Rules_Num];
        O_Act = new double[Rules_Num];
        O_Inact = new double[Rules_Num];
        
        for (int i=0; i<Rules_Num; i++) {

            min_u = 1;
            
            for (int j = 0; j< AFPRule.DescriptorsNum; j++) {
                cur_u = Rules[i].Calculate_u(j, DescValues[j]);
                if ( (cur_u!=-999) && (cur_u<min_u) )
                    min_u = cur_u;
            }
            
            u_results[i] = min_u;
            O_Act[i] = min_u * Rules[i].O_Active;
            O_Inact[i] = min_u * Rules[i].O_Inactive;
            
            Sum_u += min_u;
            Sum_Act += O_Act[i];
            Sum_Inact += O_Inact[i];
        }

        if (Sum_u!=0) {
            Final_O_Act = Sum_Act / Sum_u;
            Final_O_Inact = Sum_Inact / Sum_u;
        } else {
            Final_O_Act = 0;
            Final_O_Inact = 0;
        }
        
        
        // Calculates final response checking conditions
        
        String Response = "";
        double Max, Min;
        int ResVal = 0;
        
        if ((Final_O_Act==0)&&(Final_O_Inact==0)) {
            Response = "-";
        } else {
            if (Final_O_Act > Final_O_Inact) {
                Max = Final_O_Act;
                Min = Final_O_Inact;
                Response = "Sensitizer";
                ResVal = 2;
            } else {
                Min = Final_O_Act;
                Max = Final_O_Inact;
                Response = "NON-Sensitizer";
                ResVal = 1;
            }
            if (!((Max > Th_Accept) && ((Max-Min) > Th_Exclusion))) {
                Response = "Not Classificable";
                ResVal = -1;
            }
        }
        
        Result = new ModelAFPResult(Response, ResVal, Final_O_Act, Final_O_Inact);
        
        return Result;
    }
    
}
