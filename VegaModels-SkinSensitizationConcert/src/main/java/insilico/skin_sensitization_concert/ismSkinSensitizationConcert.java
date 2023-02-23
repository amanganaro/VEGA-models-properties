package insilico.skin_sensitization_concert;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;

public class ismSkinSensitizationConcert extends InsilicoModel {


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization_concert.xml";

    public ismSkinSensitizationConcert() throws InitFailureException {
        super(ModelData);
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
