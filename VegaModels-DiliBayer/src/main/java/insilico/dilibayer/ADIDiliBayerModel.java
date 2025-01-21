package insilico.dilibayer;

import insilico.core.ad.item.ADIndex;

public class ADIDiliBayerModel extends ADIndex {

    private double Stdev = 0;

    public ADIDiliBayerModel() {
        super("DILI model AD", "Specific AD from the Python DILI model");
    }

    public void SetStdev(double value) {
        this.Stdev = Stdev;
        SetAssessment();
    }

    @Override
    protected void SetAssessment() {
        if (this.Stdev< 0.2) {
            Assessment = "Python model standard deviation is less or equal than 0.20";
            AssessmentClass = INDEX_HIGH;
            IndexValue = 1.0;
        } else {
            Assessment = "Python model standard deviation is greater than 0.20";
            AssessmentClass = INDEX_LOW;
            IndexValue = 0.0;
        }

    }
}
