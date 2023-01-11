package insilico.melting_point_knn;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismMeltingPointKnn extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismMeltingPointKnn.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_melting_point_knn.xml";

    public ismMeltingPointKnn() throws InitFailureException {
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
