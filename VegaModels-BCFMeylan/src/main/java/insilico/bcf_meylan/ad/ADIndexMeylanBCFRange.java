package insilico.bcf_meylan.ad;

import insilico.core.ad.item.ADIndex;
import insilico.core.constant.MessagesAD;

import java.text.DecimalFormatSymbols;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class ADIndexMeylanBCFRange extends ADIndex {

    boolean IonicCompound;
    double ValueLogP;
    double ValueMW;
    
    
    public ADIndexMeylanBCFRange() {
        super(MessagesAD.RANGE_NAME, MessagesAD.RANGE_NAME_LONG);
        IonicCompound = false;
        ValueLogP = 0;
        ValueMW = 0;
    }

    
    public void setDescriptors(double LogP, double MW, boolean isIonic) {
        ValueLogP = LogP;
        ValueMW = MW;
        IonicCompound = isIonic;
        SetAssessment();
    }
    
    
    @Override
    public String GetIndexValueFormatted() {
        return IndexValue==1?"True":"False";
    }    
    
    
    @Override
    protected void SetAssessment() {
        
        DecimalFormatSymbols InternationalSymbols =
            new DecimalFormatSymbols();
        InternationalSymbols.setDecimalSeparator('.');         
        java.text.DecimalFormat df = new
        java.text.DecimalFormat("0.##", InternationalSymbols);
        
        double MW_min;
        double MW_max;
        double logP_min;
        double logP_max;

        if (IonicCompound) {
            MW_min = 68.08;
            MW_max = 991.80;
            logP_min = -6.50;
            logP_max = 11.26;
        } else {
            MW_min = 68.08;
            MW_max = 959.17;
            logP_min = -1.37;
            logP_max = 11.26;
        }
        
        if ((ValueMW < MW_min) || (ValueMW > MW_max)){
            
            Assessment = "Molecular Weight of this compound is outside the defined range [" +
                    df.format(MW_min) + "," + df.format(MW_max) + "]";
            AssessmentClass = INDEX_LOW;
            IndexValue = 0;
        
        } else if ((ValueLogP < logP_min) || (ValueLogP > logP_max)){
            
            Assessment = "LogP of this compound is outside the defined range [" +
                    df.format(logP_min) + "," + df.format(logP_max) + "]";
            AssessmentClass = INDEX_LOW;
            IndexValue = 0;

        } else {

            Assessment = "descriptors for this compound have values inside the defined range";
            AssessmentClass = INDEX_HIGH;
            IndexValue = 1;

        }             

    }
    
}
