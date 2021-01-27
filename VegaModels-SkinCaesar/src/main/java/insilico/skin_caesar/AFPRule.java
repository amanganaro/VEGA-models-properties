package insilico.skin_caesar;

/**
 * Class implementing a single rule of the AFP model
 * 
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class AFPRule {

    public final static int DescriptorsNum = 7;
    public final static int RangesNum = 4;
    
    public int Id;
    public double O_Active;
    public double O_Inactive;
    public double[][] Ranges;
    public double[][] y;
    public boolean[] Used;
    
    
    /**
     * Constructor of the class
     */
    public AFPRule() {
        
        Ranges = new double[DescriptorsNum][RangesNum];
        y = new double[DescriptorsNum][RangesNum];
        Used = new boolean[DescriptorsNum];
        for (int i=0; i<DescriptorsNum; i++)
            Used[i]=false;
        Id = 0;
        O_Active = -999;
        O_Inactive = -999;
        
    }
    

    /**
     * Sets the main parameters of the rule
     * 
     * @param Id    Rule id
     * @param O_Active  O(act) value
     * @param O_Inactive    O(inact) value
     */
    public void SetRule(int Id, double O_Active, double O_Inactive) {
        this.Id = Id;
        this.O_Active = O_Active;
        this.O_Inactive = O_Inactive;
    }

    
    /**
     * Adds coefficients of the rule for each single descriptor
     * 
     * @param DescNumber    the number of the descriptor
     * @param A     A coefficient
     * @param B     B coefficient
     * @param C     C coefficient
     * @param D     D coefficient
     * @param yA    y(A) coefficient
     * @param yB    y(B) coefficient
     * @param yC    y(C) coefficient
     * @param yD    y(D) coefficient
     */
    public void AddDescriptorRule(int DescNumber, double A, double B, double C, 
        double D, double yA, double yB, double yC, double yD) {
        Ranges[DescNumber][0] = A;
        Ranges[DescNumber][1] = B;
        Ranges[DescNumber][2] = C;
        Ranges[DescNumber][3] = D;
        y[DescNumber][0] = yA;
        y[DescNumber][1] = yB;
        y[DescNumber][2] = yC;
        y[DescNumber][3] = yD;
    }
    
    
    /**
     * Calculates the u value for a given descriptor
     * 
     * @param DescNumber    the number of descriptor
     * @param DescValue     the descriptor value
     * @return  value of u as a double
     */
    public double Calculate_u(int DescNumber, double DescValue) {
        
        if (!Used[DescNumber])
            return -999;
        
        double A = Ranges[DescNumber][0];
        double B = Ranges[DescNumber][1];
        double C = Ranges[DescNumber][2];
        double D = Ranges[DescNumber][3];
        double yA = y[DescNumber][0];
        double yB = y[DescNumber][1];
        double yC = y[DescNumber][2];
        double yD = y[DescNumber][3];
        
        if (DescValue <= A)
            return yA;
        
        if ((DescValue > A) && (DescValue < B))
            return ((yB-yA)/(B-A))*(DescValue-A) + yA;
        
        if ((DescValue >= B) && (DescValue <= C))
            return 1;
        
        if ((DescValue > C) && (DescValue < D))
            return ((yD-yC)/(D-C))*(DescValue-D) + yD;
        
        if (DescValue >= D)
            return yD;
        
        return -999;
    }
    
}
