package insilico.bcf_meylan.ad;

import insilico.core.ad.item.ADIndexADI;
import insilico.core.ad.item.iADIndex;
import insilico.core.constant.MessagesAD;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class ADIndexADIMeylanBCF extends ADIndexADI {
    
    private final double ADIValueHigh = 1.0;
    private final double ADIValueMedium = 0.85;
    private final double ADIValueLow = 0.75;
    
    
    public ADIndexADIMeylanBCF() {
        super();
        ThreshLow= 0.75;
        ThreshHigh = 0.85;
    }
    
    
    public void SetValue(double ADIValue, iADIndex AccIndex,
                         iADIndex ConcIndex, iADIndex MaxErrIndex, iADIndex LogPIndex) {
        
        if (ADIValue <= ThreshLow) {
            
            this.SetIndexValue(ADIValue);
            Assessment = MessagesAD.ADI_ASSESS_LOW;
            AssessmentClass = INDEX_LOW;
            
        } else {
            
            if ( (AccIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (ConcIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (MaxErrIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (LogPIndex.GetAssessmentClass() == INDEX_HIGH) &&
                 (ADIValue >= ThreshHigh) ) {

                Assessment = MessagesAD.ADI_ASSESS_HIGH;
                AssessmentClass = INDEX_HIGH;
                this.SetIndexValue(ADIValueHigh);
                
            } else if ( (AccIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (ConcIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (MaxErrIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (LogPIndex.GetAssessmentClass() >= INDEX_MEDIUM) &&
                 (ADIValue >= ThreshLow) ) {

                Assessment = MessagesAD.ADI_ASSESS_MEDIUM;
                AssessmentClass = INDEX_MEDIUM;
                this.SetIndexValue(ADIValueMedium);
            
            } else {

                Assessment = MessagesAD.ADI_ASSESS_LOW;
                AssessmentClass = INDEX_LOW;
                this.SetIndexValue(ADIValueLow);
                
            }
        }        
    }

}
