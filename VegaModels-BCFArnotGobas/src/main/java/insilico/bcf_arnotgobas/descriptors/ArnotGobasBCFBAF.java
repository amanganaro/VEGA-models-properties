package insilico.bcf_arnotgobas.descriptors;

/**
 * Arnot-Gobas BCFBAF equation as implemented in the EPI suite application.
 * LogKow comes from the Meylan LogP model (from EPI Suite, here implemented
 * in VEGA with some small differences) and the kM/half-life model (from EPI 
 * Suite, here implemented in VEGA with some small differences).
 * The kM value requested by the Arnot-Gobas equation is the non-log value
 * from the kM model (i.e. its main output ^ 10)
 * 
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class ArnotGobasBCFBAF {

    // Arnot-Gobas constant values
    
    static double Lbmid = 0.0685;
    static double Lblow = 0.0598;
    static double betaup = 62.7;
    static double betamid = 30.1;
    static double betalow = 16.1;
    static double Ld = 0.01;
    static double Wmid = 0.184; // middle trophic level - kg
    static double Wlow = 0.096; // lower trophic level - kg

    // Equation results

    private final double ArnotLogBCFup;
    private final double ArnotBCFup;
    private final double ArnotLogBAFup;
    private final double ArnotBAFup;

    private final double ArnotLogBCFmid;
    private final double ArnotBCFmid;
    private final double ArnotLogBAFmid;
    private final double ArnotBAFmid;

    private final double ArnotLogBCFlow;
    private final double ArnotBCFlow;
    private final double ArnotLogBAFlow;
    private final double ArnotBAFlow;
    


    public ArnotGobasBCFBAF(double UseLogKow, double HLN) {

        double Kow, phi, k1up, kDup, k2up, kEup, kGup, k1mid, kDmid, k2mid, kEmid, kGmid, k1low, kDlow, k2low, kElow, kGlow;
        double tauup, kMup, taumid, kMmid, taulow, kMlow;
        
        double Wup = 1.53; // upper trophic level - kg
        double T = 10; // deg C
        
        double Lbup = 0.107;
        
        Kow = Math.pow(10,UseLogKow);
    
        double Xpoc = 0.0000005;
        double Xdoc = 0.0000005;
        phi = 1 / (1 + (0.35 * Xpoc * Kow) + (0.08 * Xdoc * Kow));
        
        k1up = 1 / ((0.01 + 1/Kow) * Math.pow(Wup,0.4) );

        kDup = (0.02 * Math.pow(Wup,-0.15) * Math.pow(Math.E,(0.06*T))) / (0.00000005 * Kow + 2);
        k2up = k1up / (Lbup * Kow);
        kEup = 0.125 * kDup;
        kGup = 0.000502 * Math.pow(Wup, -0.2);
        kMup = 0.693/ HLN * Math.pow(Wup/0.01, -0.25);

        tauup = Math.pow((0.0065 / (((0.693/HLN) * Math.pow(0.25/0.01, -0.25)) + 0.0065)), 2);

        ArnotLogBAFup = Math.log10((1-Lbup) + (((k1up*phi) + (kDup*betaup*phi*tauup*Ld*Kow)) /
        (k2up+kEup+kGup+kMup)));

        ArnotLogBCFup = Math.log10((1-Lbup) + ((k1up*phi) / (k2up+kEup+kGup+kMup)));

        ArnotBAFup = Math.pow(10,ArnotLogBAFup);

        ArnotBCFup = Math.pow(10,ArnotLogBCFup);


        k1mid = 1 / ((0.01 + 1/Kow) * Math.pow(Wmid,0.4) );

        kDmid = (0.02 * Math.pow(Wmid,-0.15) * Math.pow(Math.E,(0.06*T))) / (0.00000005 * Kow + 2);

        k2mid = k1mid / (Lbmid * Kow);

        kEmid = 0.125 * kDmid;

        kGmid = 0.000502 * Math.pow(Wmid, -0.2);


        kMmid = 0.693/HLN * Math.pow(Wmid/0.01, -0.25);


        taumid = Math.pow((0.01 / (((0.693/HLN) * Math.pow(0.03/0.01, -0.25)) + 0.01)), 1);

        ArnotLogBAFmid = Math.log10((1-Lbmid) + (((k1mid*phi) + (kDmid*betamid*phi*taumid*Ld*Kow)) /
        (k2mid+kEmid+kGmid+kMmid)));

        ArnotLogBCFmid = Math.log10((1-Lbmid) + ((k1mid*phi) / (k2mid+kEmid+kGmid+kMmid)));

        ArnotBAFmid = Math.pow(10,ArnotLogBAFmid);

        ArnotBCFmid = Math.pow(10,ArnotLogBCFmid);


        k1low = 1 / ((0.01 + 1/Kow) * Math.pow(Wlow,0.4) );

        kDlow = (0.02 * Math.pow(Wlow,-0.15) * Math.pow(Math.E,(0.06*T))) / (0.00000005 * Kow + 2);

        k2low = k1low / (Lblow * Kow);

        kElow = 0.125 * kDlow;

        kGlow = 0.000502 * Math.pow(Wlow, -0.2);


        kMlow = 0.693/HLN * Math.pow(Wlow/0.01, -0.25);

        taulow = Math.pow((0.02 / (((0.693/HLN) * Math.pow(0.016/0.01, -0.25)) + 0.02)), 0.5);

        ArnotLogBAFlow = Math.log10((1-Lblow) + (((k1low*phi) + (kDlow*betalow*phi*taulow*Ld*Kow)) /
        (k2low+kElow+kGlow+kMlow)));

        ArnotLogBCFlow = Math.log10((1-Lblow) + ((k1low*phi) / (k2low+kElow+kGlow+kMlow)));

        ArnotBAFlow = Math.pow(10,ArnotLogBAFlow);

        ArnotBCFlow = Math.pow(10,ArnotLogBCFlow);
    }

    
    /**
     * @return the ArnotLogBCFup
     */
    public double getArnotLogBCFup() {
        return ArnotLogBCFup;
    }

    /**
     * @return the ArnotBCFup
     */
    public double getArnotBCFup() {
        return ArnotBCFup;
    }

    /**
     * @return the ArnotLogBAFup
     */
    public double getArnotLogBAFup() {
        return ArnotLogBAFup;
    }

    /**
     * @return the ArnotBAFup
     */
    public double getArnotBAFup() {
        return ArnotBAFup;
    }

    /**
     * @return the ArnotLogBCFmid
     */
    public double getArnotLogBCFmid() {
        return ArnotLogBCFmid;
    }

    /**
     * @return the ArnotBCFmid
     */
    public double getArnotBCFmid() {
        return ArnotBCFmid;
    }

    /**
     * @return the ArnotLogBAFmid
     */
    public double getArnotLogBAFmid() {
        return ArnotLogBAFmid;
    }

    /**
     * @return the ArnotBAFmid
     */
    public double getArnotBAFmid() {
        return ArnotBAFmid;
    }

    /**
     * @return the ArnotLogBCFlow
     */
    public double getArnotLogBCFlow() {
        return ArnotLogBCFlow;
    }

    /**
     * @return the ArnotBCFlow
     */
    public double getArnotBCFlow() {
        return ArnotBCFlow;
    }

    /**
     * @return the ArnotLogBAFlow
     */
    public double getArnotLogBAFlow() {
        return ArnotLogBAFlow;
    }

    /**
     * @return the ArnotBAFlow
     */
    public double getArnotBAFlow() {
        return ArnotBAFlow;
    }

}
