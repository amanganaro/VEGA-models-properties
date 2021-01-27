package insilico.skin_caesar;

/**
 * Class to store results of the AFP model
 * 
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class ModelAFPResult {

    public String Result;
    public int ResultVal;
    public double O_Active, O_Inactive;
    
    public ModelAFPResult(String Result, int ResultValue, double O_Active,
                          double O_Inactive) {
        
        this.Result = Result;
        this.ResultVal = ResultValue;
        this.O_Active = O_Active;
        this.O_Inactive = O_Inactive;
    }
    
}
