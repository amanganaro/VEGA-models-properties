package insilico.skin_sensitization_sarpy;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismSkinSensitizationSarpy extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSkinSensitizationSarpy.class);


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization_sarpy.xml";

    public ismSkinSensitizationSarpy() throws InitFailureException {
        super(ModelData);

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Activity";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        return 0;
    }

    @Override
    protected short CalculateModel() {
        return 0;
    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

    }
}
