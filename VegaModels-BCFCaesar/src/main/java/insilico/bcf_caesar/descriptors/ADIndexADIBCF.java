package insilico.bcf_caesar.descriptors;

import insilico.core.ad.item.ADIndexADI;
import insilico.core.ad.item.iADIndex;
import insilico.core.constant.MessagesAD;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class ADIndexADIBCF extends ADIndexADI {
    
    private double CIValue;
    
    
    public ADIndexADIBCF(double ADIValue, iADIndex AccIndex,
                         iADIndex ConcIndex, iADIndex MaxErrIndex) {
        
        super();
        
        if (ADIValue <= 0.75) {
            
            this.SetIndexValue(ADIValue);
            Assessment = MessagesAD.ADI_ASSESS_LOW;
            AssessmentClass = INDEX_LOW;
            
        } else {
            
            if ( (AccIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (ConcIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (MaxErrIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (ADIValue >= 0.85) ) {

                Assessment = MessagesAD.ADI_ASSESS_HIGH;
                AssessmentClass = INDEX_HIGH;
                this.SetIndexValue(1.0);
                this.CIValue = 0.76;
                return;
            }
            
            if ( (AccIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (ConcIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (MaxErrIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (ADIValue >= 0.75) ) {

                Assessment = MessagesAD.ADI_ASSESS_MEDIUM;
                AssessmentClass = INDEX_MEDIUM;
                this.SetIndexValue(0.85);
                this.CIValue = 1.12;
                return;
            }

            Assessment = MessagesAD.ADI_ASSESS_LOW;
            AssessmentClass = INDEX_LOW;
            this.SetIndexValue(0.75);
            this.CIValue = 1.27;
            
        }        
    }

    
    /**
     * @return the CIValue
     */
    public double getCIValue() {
        return CIValue;
    }
    
}
