package insilico.mitochondrial_dysfunction;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.python.CdddDescriptors;
import insilico.core.python.Communication;
import insilico.core.tools.utils.FileUtilities;
import insilico.core.tools.utils.ModelUtilities;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.JarURLConnection;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Enumeration;
import java.util.Map;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

public class MitochondrialDysfunction extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(MitochondrialDysfunction.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_mitochondrial_dysfunction.xml";

    private CdddDescriptors cdddDescriptors;

    public MitochondrialDysfunction(boolean bypassCheckCondaEnv) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData);

        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "KE2";
        this.ResultsName[1] = "AD_KE2";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        File f=File.createTempFile("output-mitochondrial-dysfunction", ".csv");
        outputTempFile = f.getAbsolutePath();

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\mitochondrial-dysfunction").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/mitochondrial-dysfunction").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            URL urlSourceEnv = MitochondrialDysfunction.class.getResource("/python/"+getCondaEnv()+".yml");
            URL urlSourceAppFile = MitochondrialDysfunction.class.getResource("/python/"+getScriptName());
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
        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        Map<String, String> Prediction = null;
        try {
            log.info("Start to execute the model");
            Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), getScriptName());

            //take the correspondent file from descriptors directory
            String descriptorFile = cdddDescriptors.getFilePathOf(CurMolecule.GetSMILES());
            Prediction=super.calculatePythonModel(pathToScriptFile, "--input "+descriptorFile,
                    " --output "+outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");

                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(ResultsName[0])));
                String[] Res = new String[ResultsSize];
                try {
                    Res[0] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(ResultsName[0])));
                } catch (Throwable ex) {
                    log.warn("Unable to find label for mitochondrial dysfunction value " + Prediction.get(ResultsName[0]));
                    Res[0] = Prediction.get(ResultsName[0]);
                }
                for(int i=1; i<ResultsSize; i++){
                    Res[i] = Prediction.get(ResultsName[i]);
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

    @Override
    public String getScriptName() {
        return "app-mitochondrial-dysfunction.py";
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
        URL urlSourceModel = getClass().getResource("/python/models-mitochondrial-dysfunction/");
        URL urlSourceDataModel = getClass().getResource("/python/data-mitochondrial-dysfunction/");

        if(urlSourceModel!=null && urlSourceEnv != null && urlSourceAppFile != null && urlSourceDataModel != null) {
            FileUtilities.copyResourcesRecursively(urlSourceModel,
                    new File(pathToExternalFolder.toString()+File.separator+"models-mitochondrial-dysfunction"));
            log.info("Models folder copied successfully");

            FileUtilities.copyResourcesRecursively(urlSourceDataModel,
                    new File(pathToExternalFolder.toString()+File.separator+"data-mitochondrial-dysfunction"));
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