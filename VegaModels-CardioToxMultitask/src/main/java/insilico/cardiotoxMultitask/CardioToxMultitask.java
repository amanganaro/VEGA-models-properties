package insilico.cardiotoxMultitask;

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
import insilico.core.tools.utils.ModelUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class CardioToxMultitask extends InsilicoModelPython {

    private static final Logger log = LoggerFactory.getLogger(CardioToxMultitask.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_cardio_tox_multitask.xml";

    private CdddDescriptors cdddDescriptors;
    private final String[] PythonResultsName;

    public CardioToxMultitask(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger) throws InitFailureException, GenericFailureException, IOException, URISyntaxException, InterruptedException {
        super(ModelData, messenger);
        isUsingCdddDescriptor=true;

        this.ResultsSize = 12;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Apical cardiotoxicity";
        this.ResultsName[1] = "Aryl hydrocarbon receptor";
        this.ResultsName[2] = "Cardiomyocyte Myocardial Injury";
        this.ResultsName[3] = "Change Action Potential";
        this.ResultsName[4] = "Change in Inotropy";
        this.ResultsName[5] = "Change In Vasoactivity";
        this.ResultsName[6] = "Endothelial injury coagulation";
        this.ResultsName[7] = "hERG channels inhibitors";
        this.ResultsName[8] = "Increase mitochondrial dysfunction";
        this.ResultsName[9] = "OxidativeStress";
        this.ResultsName[10] = "Valvular Injury Proliferation";
        this.ResultsName[11] = "ApplicabilityDomain";

        PythonResultsName = new String[ResultsSize];
        this.PythonResultsName[0] = "Apical cardiotoxicity";
        this.PythonResultsName[1] = "Aryl hydrocarbon receptor";
        this.PythonResultsName[2] = "Cardiomyocyte Myocardial Injury";
        this.PythonResultsName[3] = "Change Action Potential";
        this.PythonResultsName[4] = "Change in Inotropy";
        this.PythonResultsName[5] = "Change In Vasoactivity";
        this.PythonResultsName[6] = "Endothelial injury coagulation";
        this.PythonResultsName[7] = "hERG channels inhibitors";
        this.PythonResultsName[8] = "Increase mitochondrial dysfunction";
        this.PythonResultsName[9] = "OxidativeStress";
        this.PythonResultsName[10] = "Valvular Injury Proliferation";
        this.PythonResultsName[11] = "ApplicabilityDomain";

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        //Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),"\\AppData\\Local\\vega-models\\cardio-tox-multitask").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/cardio-tox-multitask").resolve("");
        }

        if(!bypassCheckCondaEnv) {
            boolean isEnvSet = configureCondaEnv("https://amcc.it/vega/cardio-tox-multitask.zip");
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
            Path pathToScriptFile = Paths.get(pathToExternalFolder.toString(), getScriptName());
            File f = File.createTempFile("output-cardio-tox-multitask", ".csv");
            outputTempFile = f.getAbsolutePath();
            //take the correspondent file from descriptors directory
            String descriptorFile = cdddDescriptors.getFilePathOf(CurMolecule.getInputSMILES());
            Prediction=super.calculatePythonModel(pathToScriptFile, "--input "+descriptorFile,
                    "--output " + outputTempFile);
            log.info("Finish to execute the model");

            if(Prediction != null) {
                log.info("Prediction calculated");

                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(PythonResultsName[0])));

                String[] Res = new String[ResultsSize];
                try {
                    Res[0] = this.GetTrainingSet().getClassLabel(Double.parseDouble(Prediction.get(PythonResultsName[0])));
                } catch (Throwable ex) {
                    log.warn("Unable to find label for cardio tox multitask value " + Prediction.get(PythonResultsName[0]));
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
        return "VEGA_global_V1";
    }

    @Override
    public String getScriptName() {
        return "app-cardio-tox-multitask.py";
    }

    public String getInputTempFile() {
        return inputTempFile;
    }

    @Override
    public boolean isUsingCdddDescriptor(){
        return true;
    }

}
