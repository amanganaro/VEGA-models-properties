package insilico.fathead_epa.descriptors;


/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class FatheadEPARegression {
    
    private final static double RegressionIntercept = 1.2078;
    private final static double[] RegressionCoeff = {-0.1194,
                                        -1.0817,
                                        0.0377,
                                        0.2169,
                                        0.2362,
                                        0.1303,
                                        0.1717,
                                        0.4150,
                                        0.1167,
                                        0.1974,
                                        0.1910,
                                        0.1626,
                                        0.5809,
                                        0.2251,
                                        0.6554,
                                        1.0977,
                                        -0.9756,
                                        -0.9433,
                                        0.8786,
                                        0.3271,
                                        0.6057
                                       };    
    
    
    public static double Calculate(double[] Descriptors) {
        
        double Res = RegressionIntercept;
        for (int d=0; d<RegressionCoeff.length; d++)
            Res += Descriptors[d] * RegressionCoeff[d];
        
        return Res;
    }
    
}
