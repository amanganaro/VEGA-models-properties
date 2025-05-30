package insilico.pxr_ontox;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModelPython;
import insilico.core.model.runner.iInsilicoModelRunnerMessenger;
import insilico.core.python.CdddDescriptors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class ismPxrOntox extends InsilicoModelPython {

    private static final long serialVersionUID = 1L;
    private static final String ModelData = "/data/model_pxr_ontox.xml";
    private CdddDescriptors cdddDescriptors;
    private String[] PythonResultsName;
    private static final Logger log = LoggerFactory.getLogger(ismPxrOntox.class);

    public ismPxrOntox(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger) throws InitFailureException, GenericFailureException {
        super(ModelData, messenger);
        isUsingCdddDescriptor=true;

        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "PXR prediction";
        this.ResultsName[1] = "PXR probability active";
        this.ResultsName[2] = "PXR probability non active";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "Predicted_Class";
        PythonResultsName[1] = "Probability_Active";
        PythonResultsName[2] = "Probability_NotActive";

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\pxr-ontox").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/pxr-ontox").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            boolean isEnvSet = configureCondaEnv("https://amcc.it/vega/pxr-ontox.zip");
            if(!isEnvSet) {
                throw new InitFailureException("Conda environment "+getCondaEnv()+" not set");
            }
        }

    }

    @Override
    public String getCondaEnv() {
        return "VEGA_global_V1";
    }

    @Override
    public String getScriptName() {
        return "app-pxr-ontox.py";
    }

    @Override
    public void setDescriptorGenerator(Object descriptorGenerator) {
        cdddDescriptors = (CdddDescriptors) descriptorGenerator;
    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {
            Descriptors = new double[DescriptorsSize];
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
        if(!cdddDescriptors.checkIfCdddFileIsValid(CurMolecule.getInputSMILES())) {
            return DESCRIPTORS_ERROR;
        }
        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        Map<String, String> Prediction = null;
        try {
            log.info("Start to execute the model");
            File f = File.createTempFile("output-pxr-ontox", ".csv");
            outputTempFile = f.getAbsolutePath();
            Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), getScriptName());
            String descriptorFile = cdddDescriptors.getFilePathOf(CurMolecule.getInputSMILES());

            Prediction=super.calculatePythonModel(pathToScriptFile, "--input "+descriptorFile,
                    "--output "+outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");
                double mainResult = Prediction.get(PythonResultsName[0]).equals("active") ? 1.0 : 0.0;
                CurOutput.setMainResultValue(mainResult);
                String[] Res = new String[ResultsSize];

                for(int i=0; i<PythonResultsName.length; i++){
                    try {
                        Res[i] = Format_4D.format(Prediction.get(PythonResultsName[i]));
                    } catch (Throwable ex) {
                        log.warn("Unable to find label for value " + Prediction.get(PythonResultsName[i]));
                        Res[i] = Prediction.get(PythonResultsName[i]);
                    }
                }

                CurOutput.setResults(Res);
                return MODEL_CALCULATED;
            }
            else{
                return MODEL_ERROR;
            }

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

    @Override
    public boolean isUsingCdddDescriptor(){
        return true;
    }
}
