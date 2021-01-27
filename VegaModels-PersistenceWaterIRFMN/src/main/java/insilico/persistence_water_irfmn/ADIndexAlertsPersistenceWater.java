package insilico.persistence_water_irfmn;

import insilico.core.ad.item.ADIndex;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertList;
import insilico.core.constant.MessagesAD;
import insilico.core.model.InsilicoModelOutput;

import static insilico.core.constant.InsilicoConstants.KEY_ALERT_PERS_WATER_NP;
import static insilico.core.constant.InsilicoConstants.KEY_ALERT_PERS_WATER_VP;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ADIndexAlertsPersistenceWater extends ADIndex {
    
    private static final long serialVersionUID = 1L;
    
    private final static double ThreshHigh = 1.0;
    private final static double ThreshMedium = 0.9;
    private final static double ThreshMediumNoKNN = 0.85;
    private final static double ThreshLow = 0.7;
    
    
    public ADIndexAlertsPersistenceWater() {
        super(MessagesAD.ALERT_PERSISTENCE_NAME, MessagesAD.ALERT_PERSISTENCE_NAME_LONG);
    }

    
    public void SetAlerts(AlertList FoundSAs, InsilicoModelOutput output, int UsedMols) {
        
        // check if prediction is made only on alerts
        if (UsedMols == 0) {
            SetIndexValue(ThreshMediumNoKNN);
            SetAssessment();
            return;
        }
        
        // if no alerts are available, index has a medium value
        if (FoundSAs.size() == 0) {
            SetIndexValue(ThreshMedium);
            SetAssessment();
            return;
        }
            
        // if alerts are available, check concordance
        double Prediction = output.getMainResultValue();
        int nConc=0, nNonConc=0;
        for (Alert CurAlert : FoundSAs.getSAList()) {
            if (CurAlert.getBoolProperty(KEY_ALERT_PERS_WATER_VP)) {
                if ((Prediction==3)||(Prediction==4))
                    nConc++;
                else
                    nNonConc++;
            }
            if (CurAlert.getBoolProperty(KEY_ALERT_PERS_WATER_NP)) {
                if ((Prediction==1)||(Prediction==2))
                    nConc++;
                else
                    nNonConc++;
            }
        }
        
        if (nNonConc > 0)
            SetIndexValue(ThreshLow);
        else
            SetIndexValue(ThreshHigh);
        
        SetAssessment();
    }

    @Override
    protected void SetAssessment() {
        
        if (IndexValue == ThreshHigh) {
            Assessment = MessagesAD.ALERT_PERSISTENCE_ASSESS_HIGH;
            AssessmentClass = INDEX_HIGH;
        } else if (IndexValue == ThreshMedium) {
            Assessment = MessagesAD.ALERT_PERSISTENCE_ASSESS_MEDIUM;
            AssessmentClass = INDEX_MEDIUM;
        } else if (IndexValue == ThreshMediumNoKNN) {
            Assessment = MessagesAD.ALERT_PERSISTENCE_ASSESS_MEDIUM_NO_KNN;
            AssessmentClass = INDEX_MEDIUM;
        } else {
            Assessment = MessagesAD.ALERT_PERSISTENCE_ASSESS_LOW;
            AssessmentClass = INDEX_LOW;            
        }

    }
    
}
