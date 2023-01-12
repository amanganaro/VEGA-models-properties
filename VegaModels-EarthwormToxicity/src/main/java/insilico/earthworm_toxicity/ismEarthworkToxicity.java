package insilico.earthworm_toxicity;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ismEarthworkToxicity extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismEarthworkToxicity.class);


    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_earthworm_toxicity.xml";

    public ismEarthworkToxicity() throws InitFailureException {
        super(ModelData);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "pNOEC [-log(mg/Kg)]";
        this.ResultsName[1] = "pNOEC [mg/Kg]";

        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];

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
