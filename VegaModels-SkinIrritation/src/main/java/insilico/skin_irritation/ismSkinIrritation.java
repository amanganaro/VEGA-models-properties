package insilico.skin_irritation;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismSkinIrritation extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSkinIrritation.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_irritation.xml";

    public ismSkinIrritation() throws InitFailureException {
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
