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
            case "ACHE_ONTOX":
                return "/data/model_ache_ontox.xml";
            case "AHR_ONTOX":
                return "/data/model_ahr_ontox.xml";
            case "AT1R_ONTOX":
                return "/data/model_at1r_ontox.xml";
            case "BMP_ONTOX":
                return "/data/model_bmp_ontox.xml";
            case "BSEP_ONTOX":
                return "/data/model_bsep_ontox.xml";
            case "COX1_ONTOX":
                return "/data/model_cox1_ontox.xml";
            case "CYP26_ONTOX":
                return "/data/model_cyp26_ontox.xml";
            case "FGFR1_ONTOX":
                return "/data/model_fgfr1_ontox.xml";
            case "FGFR2_ONTOX":
                return "/data/model_fgfr2_ontox.xml";
            case "FGFR3_ONTOX":
                return "/data/model_fgfr3_ontox.xml";
            case "FGFR4_ONTOX":
                return "/data/model_fgfr4_ontox.xml";
            case "GR_ONTOX":
                return "/data/model_gr_ontox.xml";
            case "HDEAC_ONTOX":
                return "/data/model_hdeac_ontox.xml";
            case "PXR_ONTOX":
                return "/data/model_pxr_ontox.xml";
            case "NMDA_ONTOX":
                return "/data/model_nmda_ontox.xml";
            case "OAT1_ONTOX":
                return "/data/model_oat1_ontox.xml";
            case "PPARA_ONTOX":
                return "/data/model_ppara_ontox.xml";
            case "PPARD_ONTOX":
                return "/data/model_ppard_ontox.xml";
            case "PPARG_ONTOX":
                return "/data/model_pparg_ontox.xml";
            case "THRA_ONTOX":
                return "/data/model_thra_ontox.xml";
            case "THRB_ONTOX":
                return "/data/model_thrb_ontox.xml";
            case "TTR_ONTOX":
                return "/data/model_ttr_ontox.xml";
            case "VGSC_ONTOX":
                return "/data/model_vgsc_ontox.xml";
            case "WNT_ONTOX":
                return "/data/model_wnt_ontox.xml";
            default:
                return "";
        }
    }

    private static String ModelFolder(String pythonModelTag) {
        switch (pythonModelTag) {
            case "ACE_ONTOX":
                return "ace-ontox_1_0_0";
            case "ACHE_ONTOX":
                return "ache-ontox_1_0_0";
            case "AHR_ONTOX":
                return "ahr-ontox_1_0_0";
            case "AT1R_ONTOX":
                return "at1r-ontox_1_0_0";
            case "BMP_ONTOX":
                return "bmp-ontox_1_0_0";
            case "BSEP_ONTOX":
                return "bsep-ontox_1_0_0";
            case "COX1_ONTOX":
                return "cox1-ontox_1_0_0";
            case "CYP26_ONTOX":
                return "cyp26-ontox_1_0_0";
            case "FGFR1_ONTOX":
                return "fgfr1-ontox_1_0_0";
            case "FGFR2_ONTOX":
                return "fgfr2-ontox_1_0_0";
            case "FGFR3_ONTOX":
                return "fgfr3-ontox_1_0_0";
            case "FGFR4_ONTOX":
                return "fgfr4-ontox_1_0_0";
            case "GR_ONTOX":
                return "gr-ontox_1_0_0";
            case "HDEAC_ONTOX":
                return "hdeac-ontox_1_0_0";
            case "PXR_ONTOX":
                return "pxr-ontox_1_0_0";
            case "NMDA_ONTOX":
                return "nmda-ontox_1_0_0";
            case "OAT1_ONTOX":
                return "nmda-oat1_1_0_0";
            case "PPARA_ONTOX":
                return "nmda-ppara_1_0_0";
            case "PPARD_ONTOX":
                return "nmda-ppard_1_0_0";
            case "PPARG_ONTOX":
                return "nmda-pparg_1_0_0";
            case "THRA_ONTOX":
                return "nmda-thra_1_0_0";
            case "THRB_ONTOX":
                return "nmda-thrb_1_0_0";
            case "TTR_ONTOX":
                return "nmda-ttr_1_0_0";
            case "VGSC_ONTOX":
                return "nmda-vgsc_1_0_0";
            case "WNT_ONTOX":
                return "nmda-wnt_1_0_0";
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

        //Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexACF().GetIndexName();

    }

    @Override
    public String getScriptName() {
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "app-ace-ontox.py";
            case "ACHE_ONTOX":
                return "app-ache-ontox.py";
            case "AHR_ONTOX":
                return "app-ahr-ontox.py";
            case "AT1R_ONTOX":
                return "app-at1r-ontox.py";
            case "BMP_ONTOX":
                return "app-bmp-ontox.py";
            case "BSEP_ONTOX":
                return "app-bsep-ontox.py";
            case "COX1_ONTOX":
                return "app-cox1-ontox.py";
            case "CYP26_ONTOX":
                return "app-cyp26-ontox.py";
            case "FGFR1_ONTOX":
                return "app-fgfr1-ontox.py";
            case "FGFR2_ONTOX":
                return "app-fgfr2-ontox.py";
            case "FGFR3_ONTOX":
                return "app-fgfr3-ontox.py";
            case "FGFR4_ONTOX":
                return "app-fgfr4-ontox.py";
            case "GR_ONTOX":
                return "app-gr-ontox.py";
            case "HDEAC_ONTOX":
                return "app-hdeac-ontox.py";
            case "PXR_ONTOX":
                return "app-pxr-ontox.py";
            case "NMDA_ONTOX":
                return "app-nmda-ontox.py";
            case "OAT1_ONTOX":
                return "app-oat1-ontox.py";
            case "PPARA_ONTOX":
                return "app-ppara-ontox.py";
            case "PPARD_ONTOX":
                return "app-ppard-ontox.py";
            case "PPARG_ONTOX":
                return "app-pparg-ontox.py";
            case "THRA_ONTOX":
                return "app-thra-ontox.py";
            case "THRB_ONTOX":
                return "app-thrb-ontox.py";
            case "TTR_ONTOX":
                return "app-ttr-ontox.py";
            case "VGSC_ONTOX":
                return "app-vgsc-ontox.py";
            case "WNT_ONTOX":
                return "app-wnt-ontox.py";
            default:
                return "";
        }
    }

    private String getReadableModelName(){
        switch(PythonModelTag){
            case "ACE_ONTOX":
                return "ACE";
            case "ACHE_ONTOX":
                return "ACHE";
            case "AHR_ONTOX":
                return "AHR";
            case "AT1R_ONTOX":
                return "AT1R";
            case "BMP_ONTOX":
                return "BMP";
            case "BSEP_ONTOX":
                return "BSEP";
            case "COX1_ONTOX":
                return "COX1";
            case "CYP26_ONTOX":
                return "CYP26";
            case "FGFR1_ONTOX":
                return "FGFR1";
            case "FGFR2_ONTOX":
                return "FGFR2";
            case "FGFR3_ONTOX":
                return "FGFR3";
            case "FGFR4_ONTOX":
                return "FGFR4";
            case "GR_ONTOX":
                return "GR";
            case "HDEAC_ONTOX":
                return "HDEAC";
            case "PXR_ONTOX":
                return "PXR";
            case "NMDA_ONTOX":
                return "NMDA";
            case "PPARA_ONTOX":
                return "PPARA";
            case "PPARD_ONTOX":
                return "PPARD";
            case "PPARG_ONTOX":
                return "PPARG";
            case "THRA_ONTOX":
                return "THRA";
            case "THRB_ONTOX":
                return "THRB";
            case "TTR_ONTOX":
                return "TTR";
            case "VGSC_ONTOX":
                return "VGSC";
            case "WNT_ONTOX":
                return "WNT";
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
