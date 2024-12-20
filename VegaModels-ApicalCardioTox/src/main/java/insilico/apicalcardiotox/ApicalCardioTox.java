package insilico.apicalcardiotox;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.model.runner.iInsilicoModelRunnerMessenger;
import insilico.core.python.CdddDescriptors;
import insilico.core.tools.utils.FileUtilities;
import insilico.core.tools.utils.ModelUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class ApicalCardioTox extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(ApicalCardioTox.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_apical_cardio_tox.xml";

    private CdddDescriptors cdddDescriptors;
    private final String[] PythonResultsName;

    public ApicalCardioTox(boolean bypassCheckCondaEnv) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Cardiotoxicity prediction";
        this.ResultsName[1] = "Python model AD assessment";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "ApicalModel";
        PythonResultsName[1] = "apical_data";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\apical-cardio-tox").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/apical-cardio-tox").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            URL urlSourceEnv = ApicalCardioTox.class.getResource("/python/"+getCondaEnv()+".yml");
            URL urlSourceAppFile = ApicalCardioTox.class.getResource("/python/"+getScriptName());
            boolean isEnvSet = configureCondaEnv(urlSourceEnv, urlSourceAppFile);
            if(!isEnvSet) {
                throw new InitFailureException("Conda environment "+getCondaEnv()+" not set");
            }
        }
    }

    public ApicalCardioTox(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData, messenger);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Cardiotoxicity prediction";
        this.ResultsName[1] = "Python model AD assessment";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "ApicalModel";
        PythonResultsName[1] = "apical_data";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\apical-cardio-tox").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/apical-cardio-tox").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            URL urlSourceEnv = ApicalCardioTox.class.getResource("/python/"+getCondaEnv()+".yml");
            URL urlSourceAppFile = ApicalCardioTox.class.getResource("/python/"+getScriptName());
            boolean isEnvSet = configureCondaEnv(urlSourceEnv, urlSourceAppFile);
            if(!isEnvSet) {
                throw new InitFailureException("Conda environment "+getCondaEnv()+" not set");
            }
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
        return DESCRIPTORS_CALCULATED;    }

    @Override
    protected short CalculateModel() {
        Map<String, String> Prediction = null;
        try {

            log.info("Start to execute the model");
            Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), getScriptName());
            File f = File.createTempFile("output-apical-cardio-tox", ".csv");
            outputTempFile = f.getAbsolutePath();
            //take the correspondent file from descriptors directory
            String descriptorFile = cdddDescriptors.getFilePathOf(CurMolecule.GetSMILES());
            Prediction=super.calculatePythonModel(pathToScriptFile, "--input "+descriptorFile,
                    "--output "+outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");

                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(PythonResultsName[0])));

                String[] Res = new String[ResultsSize];
                try {
                    Res[0] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(PythonResultsName[0])));
                } catch (Throwable ex) {
                    log.warn("Unable to find label for apical cardio tox value " + Prediction.get(PythonResultsName[0]));
                    Res[0] = Prediction.get(PythonResultsName[0]);
                }
                for(int i=1; i<ResultsSize; i++){
                    Res[i] = Prediction.get(PythonResultsName[i]);
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

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.75, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.75, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.75, 0.6);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.75, 0.6);
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == -1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        else if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
    }

    @Override
    public String getCondaEnv() {
        return "alternative";
    }

    @Override
    public String getScriptName() {
        return "app-apical-cardiotox.py";
    }

    /**
     * Add the models folder to the external path
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    @Override
    public boolean configureCondaEnv(URL urlSourceEnv, URL urlSourceAppFile) throws InterruptedException, IOException, URISyntaxException {
        boolean isSet=false;
        URL urlSourceModel = ApicalCardioTox.class.getResource("/python/models-apical-cardiotox/");
        URL urlSourceDataModel = ApicalCardioTox.class.getResource("/python/data-apical-cardiotox/");

        if(urlSourceModel!=null && urlSourceEnv != null && urlSourceAppFile != null
                && urlSourceDataModel != null) {
            FileUtilities.copyResourcesRecursively(urlSourceModel,
                    new File(pathToExternalFolder.toString()+File.separator+"models-apical-cardiotox"));
            log.info("Models folder copied successfully");

            FileUtilities.copyResourcesRecursively(urlSourceDataModel,
                    new File(pathToExternalFolder.toString()+File.separator+"data-apical-cardiotox"));
            log.info("Model data folder copied successfully");

            isSet = super.configureCondaEnv(urlSourceEnv, urlSourceAppFile);
        }
        else{
            log.error("Missing some files in setup conda {} environment", getCondaEnv());
        }
        log.info("Conda environment {} set up {}", getCondaEnv(), isSet ? "correctly": "failed");
        return isSet;
    }

    public String getInputTempFile() {
        return inputTempFile;
    }

    @Override
    public boolean isUsingCdddDescriptor(){
        return true;
    }
}
