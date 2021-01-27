package insilico.bcf_arnotgobas.descriptors.ad;

import insilico.core.ad.item.ADIndex;
import static insilico.core.ad.item.ADIndex.INDEX_HIGH;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class ADIndexLogPReliability extends ADIndex {

    private short LogPReliability;
    private boolean LogPExperimental;
    
    
    public ADIndexLogPReliability() {
        super("LogP reliability", "Reliability of logP prediction");
        LogPReliability = 0;
        LogPExperimental = false;
    }

    
    public void SetReliability(short Value, boolean isExperimental) {
        LogPReliability = Value;
        LogPExperimental = isExperimental;
        SetAssessment();
    }
    
    
    @Override
    protected void SetAssessment() {

        if ((LogPExperimental) || (LogPReliability == ADIndex.INDEX_HIGH)) {
            Assessment = "reliability of logP value used by the model is good";
            AssessmentClass = INDEX_HIGH;
            IndexValue = 1;            
        } else if (LogPReliability == ADIndex.INDEX_MEDIUM) {
            Assessment = "reliability of logP value used by the model is not optimal";
            AssessmentClass = INDEX_MEDIUM;
            IndexValue = 0.7;            
        } else if (LogPReliability == ADIndex.INDEX_LOW) {
            Assessment = "reliability of logP value used by the model is not adequate";
            AssessmentClass = INDEX_LOW;
            IndexValue = 0.0;            
        }
        
    }
    
}
