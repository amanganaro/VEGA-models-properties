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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class ismDiliBayer extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(ismDiliBayer.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_dili_bayer.xml";

    private double MainStdDev = 0;

    private CdddDescriptors cdddDescriptors;
    private String[] PythonResultsName;

    public ismDiliBayer(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData, messenger);
        isUsingCdddDescriptor=true;

        this.ResultsSize = 32;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "DILI main prediction (majority approach)";
        this.ResultsName[1] = "DILI main prediction (majority approach) st.dev";
        this.ResultsName[2] = "DILI (sensitive) prediction";
        this.ResultsName[3] = "DILI (secure) prediction";
        this.ResultsName[4] = "Prediction for assay BSEPi";
        this.ResultsName[5] = "Prediction for assay BSEPs";
        this.ResultsName[6] = "Prediction for assay PGPi";
        this.ResultsName[7] = "Prediction for assay PGPs";
        this.ResultsName[8] = "Prediction for assay MRP4i";
        this.ResultsName[9] = "Prediction for assay MRP3i";
        this.ResultsName[10] = "Prediction for assay MRP3s";
        this.ResultsName[11] = "Prediction for assay MRP2i";
        this.ResultsName[12] = "Prediction for assay MRP2s";
        this.ResultsName[13] = "Prediction for assay BCRPi";
        this.ResultsName[14] = "Prediction for assay BCRPs";
        this.ResultsName[15] = "Prediction for assay OATP1B1i";
        this.ResultsName[16] = "Prediction for assay OATP1B3i";
        this.ResultsName[17] = "Prediction for assay NRF2";
        this.ResultsName[18] = "Prediction for assay LXR";
        this.ResultsName[19] = "Prediction for assay AHR";
        this.ResultsName[20] = "Prediction for assay PPARa";
        this.ResultsName[21] = "Prediction for assay PPARg";
        this.ResultsName[22] = "Prediction for assay PXR";
        this.ResultsName[23] = "Prediction for assay FXR";
        this.ResultsName[24] = "Prediction for assay MTX_MP";
        this.ResultsName[25] = "Prediction for assay MTX_RC";
        this.ResultsName[26] = "Prediction for assay MTX_FOM";
        this.ResultsName[27] = "Prediction for assay PLD";
        this.ResultsName[28] = "Prediction for assay PLD_HTS";
        this.ResultsName[29] = "Prediction for assay HTX";
        this.ResultsName[30] = "Prediction for assay ERS";
        this.ResultsName[31] = "Prediction for assay ARE";

        PythonResultsName = new String[this.ResultsSize-1];
        PythonResultsName[0] = "DILI_majority";
        PythonResultsName[1] = "DILI_sensitive";
        PythonResultsName[2] = "DILI_secure";
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

        //Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();
        this.ADItemsName[4] = new ADIDiliBayerModel().GetIndexName();

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\dili-bayer").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/dili-bayer").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            boolean isEnvSet = configureCondaEnv("https://amcc.it/vega/dili-bayer.zip");
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
            File f = File.createTempFile("output-dili-bayer", ".csv");
            outputTempFile = f.getAbsolutePath();
            Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), getScriptName());
            String descriptorFile = cdddDescriptors.getFilePathOf(CurMolecule.getInputSMILES());

            Prediction=super.calculatePythonModel(pathToScriptFile, descriptorFile, outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");
                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(PythonResultsName[0]+"_class")));
                String[] Res = new String[ResultsSize];

                Res[0] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(PythonResultsName[0]+"_class")));
                Res[1] = Format_4D.format (Double.parseDouble(Prediction.get(PythonResultsName[0]+"_std")) );

                for(int i=1; i<PythonResultsName.length; i++){
                    try {
                        double std= Double.parseDouble(Prediction.get(PythonResultsName[i]+"_std"));
                        Res[i+1] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(PythonResultsName[i]+"_class")));
                        Res[i+1] += " - "+(std > 0.2 ? "OUT": "IN")+" AD"
                                +" (st.dev = "+ Format_4D.format(std)+")";
                    } catch (Throwable ex) {
                        log.warn("Unable to find label for value " + Prediction.get(PythonResultsName[i]+"_class"));
                        double std= Double.parseDouble(Prediction.get(PythonResultsName[i]+"_std"));
                        Res[i+1] = Prediction.get(PythonResultsName[i]+"_class");
                        Res[i+1] += " - "+(std > 0.2 ? "OUT": "IN")+" AD"
                                +" (st.dev = "+ Format_4D.format(std)+")";
                    }
                }

                MainStdDev = Double.parseDouble(Prediction.get(PythonResultsName[0]+"_std"));

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

        // Specific python AD
        ADIDiliBayerModel ADpython = new ADIDiliBayerModel();
        ADpython.SetStdev(MainStdDev);
        this.CurOutput.addADIndex(ADpython);

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;
        if (ADpython.GetIndexValue() == 0.0)
            ADIValue = 0.0;

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
        return "VEGA_liver_mtnn";
    }

    @Override
    public String getScriptName() {
        return "app-dili-bayer.py";
    }

    public String getInputTempFile() {
        return inputTempFile;
    }

    @Override
    public boolean isUsingCdddDescriptor(){
        return true;
    }
}
