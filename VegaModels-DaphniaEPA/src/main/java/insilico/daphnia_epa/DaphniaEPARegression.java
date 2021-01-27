package insilico.daphnia_epa;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class DaphniaEPARegression {
    
    private final static double RegressionIntercept = -6.5594;
    private final static double[] RegressionCoeff = {1.2157,
                                        0.1341,
                                        0.6974,
                                        -1.3213,
                                        0.8605,
                                        1.4685,
                                        -0.9197,
                                        0.2238,
                                        1.4502,
                                        2.4060,
                                        1.9085,
                                        -2.4036,
                                        -0.3463,
                                        0.0255,
                                        -1.4215,
                                        -0.7185,
                                        -1.0232,
                                        -1.5228
                                       };
    
    
    public static double Calculate(double[] Descriptors) {
        
        double Res = RegressionIntercept;
        for (int d=0; d<RegressionCoeff.length; d++)
            Res += Descriptors[d] * RegressionCoeff[d];
        
        return Res;
    }
    
}
