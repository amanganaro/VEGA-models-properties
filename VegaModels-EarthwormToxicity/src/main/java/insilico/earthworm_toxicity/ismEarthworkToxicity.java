package insilico.earthworm_toxicity;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.earthworm_toxicity.descriptors.EmbeddedDescriptors;
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
        this.ResultsName[1] = "NOEC [mg/Kg]";

        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "B06[N−O]";
        DescriptorsNames[1] = "SaasC";
        DescriptorsNames[2] = "B06[N−F]";
        DescriptorsNames[3] = "GATS1m";
        DescriptorsNames[4] = "VE1signB(e)";
        DescriptorsNames[5] = "B03[N−Cl]";
        DescriptorsNames[6] = "B08[N−Cl]";
        DescriptorsNames[7] = "S−107";


    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.B06_N_O;
            Descriptors[1] = embeddedDescriptors.SaasC;
            Descriptors[2] = embeddedDescriptors.B06_N_F;
            Descriptors[3] = embeddedDescriptors.GATS1m;
            Descriptors[4] = embeddedDescriptors.VE1signB_e;
            Descriptors[5] = embeddedDescriptors.B03_N_Cl;
            Descriptors[6] = embeddedDescriptors.B08_N_Cl;
            Descriptors[7] = embeddedDescriptors.S_107;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        double Prediction;
        try {

            Prediction = 1.062
                        - 1.389 * Descriptors[0]
                        - 0.156 * Descriptors[1]
                        + 1.29 * Descriptors[2]
                        - 1.9 * Descriptors[3]
                        - 2.425 * Descriptors[4]
                        - 0.821 * Descriptors[5]
                        + 0.225 * Descriptors[6]
                        + 1.233 * Descriptors[7];

            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = String.valueOf(Format_3D.format(Prediction));
            Res[1] = String.valueOf(Format_3D.format(Math.pow(10,-Prediction)));
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex){
            return MODEL_ERROR;
        }

    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

    }
}
