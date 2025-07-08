package insilico.ontox_assay;

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
import org.apache.commons.lang3.math.NumberUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class ismOntoxAssay extends InsilicoModelPython {

    private static final long serialVersionUID = 1L;
    private CdddDescriptors cdddDescriptors;
    private String[] PythonResultsName;
    private static final Logger log = LoggerFactory.getLogger(ismOntoxAssay.class);
    private String PythonModelTag = "";

    private static String ModelData(String pythonModelTag) {
        switch (pythonModelTag) {
            case "ACE_ONTOX":
                return "/data/model_ace_ontox.xml";
            case "PXR_ONTOX":
                return "/data/model_pxr_ontox.xml";
            case "NMDA_ONTOX":
                return "/data/model_nmda_ontox.xml";
            default:
                return "";
        }
    }

    private static String ModelFolder(String pythonModelTag) {
        switch (pythonModelTag) {
            case "ACE_ONTOX":
                return "ace-ontox_1_0_0";
            case "PXR_ONTOX":
                return "pxr-ontox_1_0_0";
            case "NMDA_ONTOX":
                return "nmda-ontox_1_0_0";
            default:
                return "";
        }
    }

    public ismOntoxAssay(boolean bypassCheckCondaEnv, iInsilicoModelRunnerMessenger messenger, String pythonModelTag) throws InitFailureException, GenericFailureException {
        super(ModelData(pythonModelTag), messenger, ModelFolder(pythonModelTag), "GLOBAL", bypassCheckCondaEnv);

        PythonModelTag = pythonModelTag;
        isUsingCdddDescriptor=true;

        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = getReadableModelName()+" prediction";
        this.ResultsName[1] = getReadableModelName()+" probability active";
        this.ResultsName[2] = getReadableModelName()+" probability non active";

        PythonResultsName = new String[this.ResultsSize];
        PythonResultsName[0] = "Predicted_Class";
        PythonResultsName[1] = "Probability_Active";
        PythonResultsName[2] = "Probability_NotActive";

    }

    @Override
    public String getScriptName() {
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "app-ace-ontox.py";
            case "PXR_ONTOX":
                return "app-pxr-ontox.py";
            case "NMDA_ONTOX":
                return "app-nmda-ontox.py";
            default:
                return "";
        }
    }

    private String getReadableModelName(){
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "ACE";
            case "PXR_ONTOX":
                return "PXR";
            case "NMDA_ONTOX":
                return "NMDA";
            default:
                return "";
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
            File f = File.createTempFile("output-ontox-assay", ".csv");
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
                        if(NumberUtils.isCreatable(Prediction.get(PythonResultsName[i]))){
                            Res[i] = Format_4D.format(Double.parseDouble(Prediction.get(PythonResultsName[i])));
                        }
                        else{
                            double pred_number =  Prediction.get(PythonResultsName[i]).equals("active") ? 1.0 : 0.0;
                            Res[i] = this.GetTrainingSet().getClassLabel(pred_number);
                        }

                    } catch (Throwable ex) {
                        log.warn("Unable to find label for value {}", Prediction.get(PythonResultsName[i]));
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
    public boolean isUsingCdddDescriptor(){
        return true;
    }
}
