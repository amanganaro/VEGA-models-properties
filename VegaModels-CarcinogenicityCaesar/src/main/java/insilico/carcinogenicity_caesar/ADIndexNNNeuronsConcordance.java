package insilico.carcinogenicity_caesar;

import insilico.core.ad.item.ADIndex;

/**
 *
 * @author amanganaro
 */
public class ADIndexNNNeuronsConcordance extends ADIndex {

    private static final long serialVersionUID = 1L;
    
    
    public ADIndexNNNeuronsConcordance() {
        super("Neurons concordance", "Neural map neurons concordance");
        this.DecimalDigits = 2;
    }
    
    
    public void SetNeuronData(int Population, double Concordance) {
        
        if (Population == 0) {
            this.IndexValue = 0.5;
        } else  {
            if (Concordance > 0.65) {
                this.IndexValue = 1;
            } else {
                this.IndexValue = 0.75;
            }
        }
        
        SetAssessment();
    }
    
    
    @Override
    protected void SetAssessment() {
        if (this.IndexValue == 1) {
            Assessment = "predicted value agrees with experimental values of "
                        + "training set compounds laying in the same neuron";
            AssessmentClass = INDEX_HIGH;
        } else if (this.IndexValue == 0.75) {
            Assessment = "predicted value disagrees with experimental values of "
                        + "training set compounds laying in the same neuron";
            AssessmentClass = INDEX_MEDIUM;
        } else if (this.IndexValue == 0.5) {
            Assessment = "predicted substance falls into a neuron that is populated by "
                        + " no compounds of the training set";
            AssessmentClass = INDEX_LOW;
        } 
    }
    
}
