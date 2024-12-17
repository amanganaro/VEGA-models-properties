package insilico.dilibayer;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.python.CdddDescriptors;
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

public class ismDiliBayer extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(ismDiliBayer.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_dili_bayer.xml";

    private CdddDescriptors cdddDescriptors;

    public ismDiliBayer(boolean bypassCheckCondaEnv) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData);

        this.ResultsSize = 31;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "DILI_secure";
        this.ResultsName[1] = "DILI_sensitive";
        this.ResultsName[2] = "DILI_majority";
        this.ResultsName[3] = "BSEPi";
        this.ResultsName[4] = "BSEPs";
        this.ResultsName[5] = "PGPi";
        this.ResultsName[6] = "PGPs";
        this.ResultsName[7] = "MRP4i";
        this.ResultsName[8] = "MRP3i";
        this.ResultsName[9] = "MRP3s";
        this.ResultsName[10] = "MRP2i";
        this.ResultsName[11] = "MRP2s";
        this.ResultsName[12] = "BCRPi";
        this.ResultsName[13] = "BCRPs";
        this.ResultsName[14] = "OATP1B1i";
        this.ResultsName[15] = "OATP1B3i";
        this.ResultsName[16] = "NRF2";
        this.ResultsName[17] = "LXR";
        this.ResultsName[18] = "AHR";
        this.ResultsName[19] = "PPARa";
        this.ResultsName[20] = "PPARg";
        this.ResultsName[21] = "PXR";
        this.ResultsName[22] = "FXR";
        this.ResultsName[23] = "MTX_MP";
        this.ResultsName[24] = "MTX_RC";
        this.ResultsName[25] = "MTX_FOM";
        this.ResultsName[26] = "PLD";
        this.ResultsName[27] = "PLD_HTS";
        this.ResultsName[28] = "HTX";
        this.ResultsName[29] = "ERS";
        this.ResultsName[30] = "ARE";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        File f = File.createTempFile("output-dili-bayer", ".csv");
        outputTempFile = f.getAbsolutePath();

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\dili-bayer").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/dili-bayer").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            URL urlSourceEnv = ismDiliBayer.class.getResource("/python/"+getCondaEnv()+".yml");
            URL urlSourceAppFile = ismDiliBayer.class.getResource("/python/"+getScriptName());
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

            Prediction=super.calculatePythonModel(pathToScriptFile, descriptorFile, outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");
                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(ResultsName[0]+"_class")));
                String[] Res = new String[ResultsSize];

                for(int i=0; i<ResultsSize; i++){
                    try {
                        Res[i] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(ResultsName[i]+"_class")));
                    } catch (Throwable ex) {
                        log.warn("Unable to find label for carcinogenicity value " + Prediction.get(ResultsName[i]+"_class"));
                        Res[i] = Prediction.get(ResultsName[i]+"_class");
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
        return "liver-mtnn";
    }

    @Override
    public String getScriptName() {
        return "app-dili-bayer.py";
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
        URL urlSourceModel = getClass().getResource("/python/models_dili_bayer/");

        if(urlSourceModel!=null && urlSourceEnv != null && urlSourceAppFile != null){
            FileUtilities.copyResourcesRecursively(urlSourceModel,
                    new File(pathToExternalFolder.toString()+File.separator+"models_dili_bayer"));
            log.info("Models folder copied successfully");

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
