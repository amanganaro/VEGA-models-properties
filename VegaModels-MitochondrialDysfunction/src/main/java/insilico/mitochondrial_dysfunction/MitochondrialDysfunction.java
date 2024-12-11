package insilico.mitochondrial_dysfunction;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.tools.utils.FileUtilities;
import insilico.core.tools.utils.ModelUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class MitochondrialDysfunction extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(MitochondrialDysfunction.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_mitochondrial_dysfunction.xml";

    protected final boolean CHECK_SETUP = true;

    protected String descriptorsTempFile = "";

    public MitochondrialDysfunction() throws InitFailureException, GenericFailureException, IOException {
        super(ModelData);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "KE2";
        this.ResultsName[1] = "AD_KE2";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        File f = File.createTempFile("input-mitochondrial-dysfunction", ".csv");
        inputTempFile = f.getAbsolutePath();
        f=File.createTempFile("output-mitochondrial-dysfunction", ".csv");
        outputTempFile = f.getAbsolutePath();
        f=File.createTempFile("descriptors-mitochondrial-dysfunction", ".csv");
        descriptorsTempFile = f.getAbsolutePath();

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\mitochondrial-dysfunction\\python").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/mitochondrial-dysfunction/python").resolve("");
        }
    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        log.info("enter in the calculate descriptors method");
        try {
            prepareInputData();
            Descriptors cdddDescriptors = new Descriptors(this);
            boolean result=cdddDescriptors.calculateDescriptors();

            if(!result){
                log.info("Descriptors calculation failed");
                return DESCRIPTORS_ERROR;
            }

            log.info("Descriptors calculated correctly");
            Descriptors = new double[DescriptorsSize];

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        log.info("enter in the calculate model method");
        Map<String, String> Prediction;
        try {
            boolean isEnvSet = CHECK_SETUP ? configureCondaEnv() : true;
            if(isEnvSet){
                log.info("Start to execute the model");
                Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), "app.py");
                Prediction=super.calculatePythonModel(pathToScriptFile, "--input "+descriptorsTempFile,
                        "--output "+outputTempFile);
                log.info("Finish to execute the model");
            }
            else{
                Prediction = null;
            }

            if(Prediction != null) {
                log.info("Prediction calculated");

                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(ResultsName[0])));

                String[] Res = new String[ResultsSize];
                for(int i=0; i<ResultsSize; i++){
                    Res[i] = Prediction.get(ResultsName[i]);
                }

                CurOutput.setResults(Res);

                File file = new File(descriptorsTempFile);
                file.delete();

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

        // questo lo vediamo quando troviamo il dataset da inserire
        return InsilicoModel.AD_ERROR;

//        // Calculates various AD indices
//        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
//        adq.setMoleculesForIndexSize(2);
//        if (!adq.Calculate(CurMolecule, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets threshold for AD indices
//        try {
//            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
//            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.6);
//            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.6);
//            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.6);
//        } catch (Throwable e) {
//            return InsilicoModel.AD_ERROR;
//        }
//
//        // Sets Range check
//        ADCheckDescriptorRange adrc = new ADCheckDescriptorRange();
//        if (!adrc.Calculate(TS, insilico.dilibayer.Descriptors, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets ACF check
//        ADCheckACF adacf = new ADCheckACF(TS);
//        if (!adacf.Calculate(CurMolecule, CurOutput))
//            return InsilicoModel.AD_ERROR;
//
//        // Sets final AD index
//        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
//        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
//        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;
//
//        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
//        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
//                CurOutput.getADIndex(ADIndexConcordance.class),
//                CurOutput.getADIndex(ADIndexMaxError.class));
//        CurOutput.setADI(ADI);
//
//        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }

    @Override
    public String getCondaEnv() {
        return "alternative";
    }

    /**
     * Add the models folder to the external path
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    public boolean configureCondaEnv() throws IOException, InterruptedException, URISyntaxException {

        log.info("enter in the configure conda env method");

        boolean isSet=false;
        URL urlSourceEnv = getClass().getResource("/python/"+getCondaEnv()+".yml");
        URL urlSourceAppFile = getClass().getResource("/python/app.py");
        URL urlSourceModel = getClass().getResource("/python/models/");
        URL urlSourceDataModel = getClass().getResource("/python/data/");

        if(urlSourceModel!=null && urlSourceEnv != null && urlSourceAppFile != null
                && urlSourceDataModel != null) {
            FileUtilities.copyExternalData(Paths.get(urlSourceModel.toURI()).toString(),
                    (pathToExternalFolder.toString()+File.separator+"models"));
            log.info("Models folder copied successfully");

            FileUtilities.copyExternalData(Paths.get(urlSourceDataModel.toURI()).toString(),
                    (pathToExternalFolder.toString()+File.separator+"data"));
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

    public String getDescriptorsTempFile() {
        return descriptorsTempFile;
    }
}
