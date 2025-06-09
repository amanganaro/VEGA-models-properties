package insilico.pmml_ontox;

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

public class ismPmmlOntox extends InsilicoModelPython {

    private static final long serialVersionUID = 1L;
    //private static String ModelData = "/data/model_ace_ontox.xml";
    private CdddDescriptors cdddDescriptors;
    private String[] PythonResultsName;
    private static final Logger log = LoggerFactory.getLogger(ismPmmlOntox.class);
    private String PythonModelTag = "";


    private static String ModelData(String pythonModelTag) {
        switch (pythonModelTag) {
            case "ACE_ONTOX":
                return "/data/model_ace_ontox.xml";
            case "TEST_ONTOX":
                return "/data/model_test_ontox.xml";
            default:
                return "/data/model_ace_ontox.xml";
        }
    }

    public ismPmmlOntox(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger, String pythonModelTag) throws InitFailureException, GenericFailureException {
        super(ModelData(pythonModelTag), messenger);

        PythonModelTag = pythonModelTag;
        isUsingCdddDescriptor=true;

        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = getReadableModelName()+" prediction";
        this.ResultsName[1] = getReadableModelName()+" probability";
        this.ResultsName[2] = getReadableModelName()+" probability active";
        this.ResultsName[3] = getReadableModelName()+" probability non active";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "predicted_Category";
        PythonResultsName[1] = "probability";
        PythonResultsName[2] = "probability_active";
        PythonResultsName[3] = "probability_not active";

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\pmml-ontox").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/pmml-ontox").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            boolean isEnvSet = configureCondaEnv("https://amcc.it/vega/pmml-ontox.zip");
            if(!isEnvSet) {
                throw new InitFailureException("Conda environment "+getCondaEnv()+" not set");
            }
        }
    }


    @Override
    public String getCondaEnv() {
        return "VEGA_global_V2";
    }

    @Override
    public String getScriptName() {
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "app-ace-ontox.py";
            case "TEST_ONTOX":
                return "app-test-ontox.py";
            default:
                return "app-ace-ontox.py";
        }
    }

    private String getReadableModelName(){
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "ACE";
            case "TEST_ONTOX":
                return "Test";
            default:
                return "Ace";
        }
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
            File f = File.createTempFile("output-pmml-ontox", ".csv");
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
