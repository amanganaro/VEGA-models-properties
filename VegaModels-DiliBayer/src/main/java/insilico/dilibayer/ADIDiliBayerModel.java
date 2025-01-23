package insilico.dilibayer;

import insilico.core.ad.item.ADIndex;

public class ADIDiliBayerModel extends ADIndex {

    private double Stdev = 0;

    public ADIDiliBayerModel() {
        super("Standard deviation check", "Check of the standard deviation of the ensemble model");
    }

    public void SetStdev(double value) {
        this.Stdev = value;
        SetAssessment();
    }

    @Override
    protected void SetAssessment() {
        if (this.Stdev< 0.2) {
            Assessment = "Standard deviation of ensemble model probability is less or equal than 0.20";
            AssessmentClass = INDEX_HIGH;
            IndexValue = 1.0;
        } else {
            Assessment = "Standard deviation of ensemble model probability is greater than 0.20";
            AssessmentClass = INDEX_LOW;
            IndexValue = 0.0;
        }

    }
}
