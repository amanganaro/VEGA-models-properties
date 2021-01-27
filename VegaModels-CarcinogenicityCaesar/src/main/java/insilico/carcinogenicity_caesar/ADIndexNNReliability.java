package insilico.carcinogenicity_caesar;

import insilico.core.ad.item.ADIndex;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ADIndexNNReliability extends ADIndex {
    
    private static final long serialVersionUID = 1L;
    
    private final static double Threshold = 0.1;
    private double Delta;
    
    public ADIndexNNReliability() {
        super("Pos/Non-Pos difference", "Model class assignment reliability");
        this.Delta = 0;
    }

    
    public void SetDelta(double Delta) {
        this.Delta = Delta;
        SetAssessment();
    }

    @Override
    protected void SetAssessment() {
        this.IndexValue = this.Delta;
        if (Delta > Threshold) {
            Assessment = "model class assignment is well defined";
            AssessmentClass = INDEX_HIGH;
        } else {
            Assessment = "model class assignment is uncertain";
            AssessmentClass = INDEX_LOW;
        } 
    }

    
}
