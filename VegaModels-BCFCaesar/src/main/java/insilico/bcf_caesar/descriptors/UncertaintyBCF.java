package insilico.bcf_caesar.descriptors;

import insilico.core.ad.reasoning.Uncertainty;
import insilico.core.ad.reasoning.UncertaintyClassBar;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class UncertaintyBCF extends Uncertainty {
    
    
    public UncertaintyBCF() {
        super();
    }
    
    
    public void Evaluate(double Prediction, double ADI) {
        
        DecimalFormatSymbols InternationalSymbols = new DecimalFormatSymbols();
        InternationalSymbols.setDecimalSeparator('.');         
        DecimalFormat df = new DecimalFormat("0.##", InternationalSymbols);
        
        
        
        //// 3.3 Threshold
        
        UncertaintyClassBar Bar_3_3 = new UncertaintyClassBar();
        
        Bar_3_3.setClassName("Threshold 3.3 (bioaccumulative)");
        Bar_3_3.setAxisName("logBCF");
        Bar_3_3.setXValue(Prediction);
        Bar_3_3.setXValueAsString(df.format(Prediction));
        Bar_3_3.getThresholdsMarks().add("3.3 (B)");
        Bar_3_3.getThresholds().add(3.3);

        if (ADI > 0.85) {
            Bar_3_3.setXValueInterval(0.5);
        } else if (ADI > 0.75) {
            Bar_3_3.setXValueInterval(0.7);
        } else {
            Bar_3_3.setXValueInterval(0.0);
        }
        
        if (Bar_3_3.getXValueInterval() > 0) {
            double ClassVal = Prediction + Bar_3_3.getXValueInterval();
            String Classification = "";
            if (ClassVal < 3.3) 
                Classification = "can be safely classified as not bioaccumulative";
            else if (ClassVal > 3.3) {
                if (Prediction <= 3.3)
                    Classification = "can not be safely classified as not bioaccumulative";
                else
                    Classification = "can be classified as bioaccumulative";
            }

            Bar_3_3.setClassDescription("For the threshold logBCF = 3.3, the current compound can be associated (due to its Applicability Domain index value) " +
                    "to a conservative interval of " + Bar_3_3.getXValueInterval() + " log units. \nOn this basis, the compound " +
                    Classification + ".");
        } else {
            Bar_3_3.setClassDescription("For the threshold logBCF = 3.3, the current compound can not be associated (due to its Applicability Domain index value) " +
                    "to any conservative interval. \nNo safe classification can be done.");
        }

        this.Bars.add(Bar_3_3);

        
        
        //// 3.7 Threshold
        
        UncertaintyClassBar Bar_3_7 = new UncertaintyClassBar();
        
        Bar_3_7.setClassName("Threshold 3.7 (very bioaccumulative)");
        Bar_3_7.setAxisName("logBCF");
        Bar_3_7.setXValue(Prediction);
        Bar_3_7.setXValueAsString(df.format(Prediction));
        Bar_3_7.getThresholdsMarks().add("3.7 (vB)");
        Bar_3_7.getThresholds().add(3.7);

        if (ADI > 0.85) {
            Bar_3_7.setXValueInterval(0.5);
        } else if (ADI > 0.75) {
            Bar_3_7.setXValueInterval(0.5);
        } else {
            Bar_3_7.setXValueInterval(0.0);
        }
        
        if (Bar_3_7.getXValueInterval() > 0) {
            double ClassVal = Prediction + Bar_3_7.getXValueInterval();
            String Classification = "";
            if (ClassVal < 3.7) 
                Classification = "can be safely classified as not very bioaccumulative";
            else if (ClassVal > 3.7) {
                if (Prediction <= 3.7)
                    Classification = "can not be safely classified as not very bioaccumulative";
                else
                    Classification = "can be classified as very bioaccumulative";
            }

            Bar_3_7.setClassDescription("For the threshold logBCF = 3.7, the current compound can be associated (due to its Applicability Domain index value) " +
                    "to a conservative interval of " + Bar_3_7.getXValueInterval() + " log units. \nOn this basis, the compound " +
                    Classification + ".");
        } else {
            Bar_3_7.setClassDescription("For the threshold logBCF = 3.7, the current compound can not be associated (due to its Applicability Domain index value) " +
                    "to any conservative interval. \nNo safe classification can be done.");
        }

        this.Bars.add(Bar_3_7);        
    }    
    
}
