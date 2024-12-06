import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.InsilicoModelPython;
import insilico.core.tools.utils.FileUtilities;
import insilico.core.tools.utils.ModelUtilities;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

public class ismDiliBayer extends InsilicoModelPython {

    private static final Logger log = LogManager.getLogger(ismDiliBayer.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_dili_bayer.xml";

    public ismDiliBayer() throws InitFailureException {
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

    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {
            //TEMPPPP to put into pythonSilicoModel
            //save into input.csv the smiles to be processed
            FileWriter myWriter = new FileWriter("input.csv");
            myWriter.write("smiles\r\n"+CurMolecule.GetSMILES());
            myWriter.close();

            Descriptors cdddDescriptors = new Descriptors(CurMolecule.GetSMILES());
            boolean result=cdddDescriptors.calculateDescriptors();

            if(!result){
                return DESCRIPTORS_ERROR;
            }

            Descriptors = new double[DescriptorsSize];

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        Map<String, String> Prediction;
        try {

            boolean isEnvSet = configureCondaEnv();
            if(isEnvSet){
                Prediction = super.calculatePythonModel("DILI_secure_mean");
            }
            else{
                Prediction = null;
            }


            if(Prediction != null) {

                CurOutput.setMainResultValue(Double.parseDouble(Prediction.get(ResultsName[0])));

                String[] Res = new String[ResultsSize];
                for(int i=0; i<ResultsSize; i++){
                    Res[i] = Prediction.get(ResultsName[i]+"_class");
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
//        if (!adrc.Calculate(TS, Descriptors, CurOutput))
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

    public boolean configureCondaEnv() throws IOException, InterruptedException {
        boolean isSet=false;

        try {
            URL resource = ismDiliBayer.class.getResource("python" + File.separator + getCondaEnv()+".yml");
            if(resource != null){
                Path pathToEnvFile= Paths.get(resource.toURI()).toAbsolutePath();
                isSet = super.configureCondaEnv(pathToEnvFile);
            }
            else{
                log.error("Cannot find file {}", getCondaEnv()+".yml");
                return false;
            }
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        }

        return isSet;
    }
}
