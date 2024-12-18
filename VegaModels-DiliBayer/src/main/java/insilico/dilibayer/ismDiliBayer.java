package insilico.dilibayer;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
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
    private final String[] PythonResultsName;

    public ismDiliBayer(boolean bypassCheckCondaEnv) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData);

        this.ResultsSize = 31;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "DILI (secure) prediction";
        this.ResultsName[1] = "DILI (sensitive) prediction";
        this.ResultsName[2] = "DILI (majority) prediction";
        this.ResultsName[3] = "Prediction for essay BSEPi";
        this.ResultsName[4] = "Prediction for essay BSEPs";
        this.ResultsName[5] = "Prediction for essay PGPi";
        this.ResultsName[6] = "Prediction for essay PGPs";
        this.ResultsName[7] = "Prediction for essay MRP4i";
        this.ResultsName[8] = "Prediction for essay MRP3i";
        this.ResultsName[9] = "Prediction for essay MRP3s";
        this.ResultsName[10] = "Prediction for essay MRP2i";
        this.ResultsName[11] = "Prediction for essay MRP2s";
        this.ResultsName[12] = "Prediction for essay BCRPi";
        this.ResultsName[13] = "Prediction for essay BCRPs";
        this.ResultsName[14] = "Prediction for essay OATP1B1i";
        this.ResultsName[15] = "Prediction for essay OATP1B3i";
        this.ResultsName[16] = "Prediction for essay NRF2";
        this.ResultsName[17] = "Prediction for essay LXR";
        this.ResultsName[18] = "Prediction for essay AHR";
        this.ResultsName[19] = "Prediction for essay PPARa";
        this.ResultsName[20] = "Prediction for essay PPARg";
        this.ResultsName[21] = "Prediction for essay PXR";
        this.ResultsName[22] = "Prediction for essay FXR";
        this.ResultsName[23] = "Prediction for essay MTX_MP";
        this.ResultsName[24] = "Prediction for essay MTX_RC";
        this.ResultsName[25] = "Prediction for essay MTX_FOM";
        this.ResultsName[26] = "Prediction for essay PLD";
        this.ResultsName[27] = "Prediction for essay PLD_HTS";
        this.ResultsName[28] = "Prediction for essay HTX";
        this.ResultsName[29] = "Prediction for essay ERS";
        this.ResultsName[30] = "Prediction for essay ARE";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "DILI_secure";
        PythonResultsName[1] = "DILI_sensitive";
        PythonResultsName[2] = "DILI_majority";
        PythonResultsName[3] = "BSEPi";
        PythonResultsName[4] = "BSEPs";
        PythonResultsName[5] = "PGPi";
        PythonResultsName[6] = "PGPs";
        PythonResultsName[7] = "MRP4i";
        PythonResultsName[8] = "MRP3i";
        PythonResultsName[9] = "MRP3s";
        PythonResultsName[10] = "MRP2i";
        PythonResultsName[11] = "MRP2s";
        PythonResultsName[12] = "BCRPi";
        PythonResultsName[13] = "BCRPs";
        PythonResultsName[14] = "OATP1B1i";
        PythonResultsName[15] = "OATP1B3i";
        PythonResultsName[16] = "NRF2";
        PythonResultsName[17] = "LXR";
        PythonResultsName[18] = "AHR";
        PythonResultsName[19] = "PPARa";
        PythonResultsName[20] = "PPARg";
        PythonResultsName[21] = "PXR";
        PythonResultsName[22] = "FXR";
        PythonResultsName[23] = "MTX_MP";
        PythonResultsName[24] = "MTX_RC";
        PythonResultsName[25] = "MTX_FOM";
        PythonResultsName[26] = "PLD";
        PythonResultsName[27] = "PLD_HTS";
        PythonResultsName[28] = "HTX";
        PythonResultsName[29] = "ERS";
        PythonResultsName[30] = "ARE";


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
                        Res[i] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(PythonResultsName[i]+"_class")));
                    } catch (Throwable ex) {
                        log.warn("Unable to find label for value " + Prediction.get(PythonResultsName[i]+"_class"));
                        Res[i] = Prediction.get(PythonResultsName[i]+"_class");
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
