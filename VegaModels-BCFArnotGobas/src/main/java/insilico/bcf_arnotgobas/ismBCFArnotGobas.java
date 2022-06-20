package insilico.bcf_arnotgobas;

import insilico.bcf_arnotgobas.descriptors.KmFactor;
import insilico.bcf_arnotgobas.descriptors.MeylanLogP;
import insilico.bcf_arnotgobas.descriptors.ad.ADIndexADIMeylanBCF;
import insilico.bcf_arnotgobas.descriptors.ad.ADIndexLogPReliability;
import insilico.bcf_arnotgobas.descriptors.ArnotGobasBCFBAF;
import insilico.bcf_arnotgobas.descriptors.DescriptorAnalysisBCFMeylan;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import insilico.km_arnot.ismKmArnot;
import insilico.meylanlogp.ismLogPMeylan;


import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismBCFArnotGobas extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_bcf_arnotgobas.xml";

    private final ismLogPMeylan LogPModel;
    private final ismKmArnot kMModel;

    private InsilicoModelOutput LogPResult;
    private double LogPValueForModel;
    private InsilicoModelOutput kMResult;
    private double kMValueForModel;


    public ismBCFArnotGobas()
            throws InitFailureException {
        super(ModelData);

        // Builds model objects
        LogPModel = new ismLogPMeylan();
        kMModel = new ismKmArnot();

        // Define no. of descriptors
        this.DescriptorsSize = 2;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "LogP";
        this.DescriptorsNames[1] = "kM";

        // Defines results
        this.ResultsSize = 10;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted BCF (up) [log(L/kg)]";
        this.ResultsName[1] = "Predicted BCF (up) [L/kg]";
        this.ResultsName[2] = "Predicted BCF (low) [log(L/kg)]";
        this.ResultsName[3] = "Predicted BCF (low) [L/kg]";
        this.ResultsName[4] = "Predicted BCF (mid) [log(L/kg)]";
        this.ResultsName[5] = "Predicted BCF (mid) [L/kg]";
        this.ResultsName[6] = "Predicted LogP (Meylan/Kowwin)";
        this.ResultsName[7] = "Predicted LogP reliability";
        this.ResultsName[8] = "Predicted kM (Meylan)";
        this.ResultsName[9] = "Predicted kM reliability";

        // Define AD items
        this.ADItemsName = new String[7];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexLogPReliability().GetIndexName();
        this.ADItemsName[5] = new ADIndexACF().GetIndexName();
        // aggiungere reliability per kM?

    }


    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {

        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        DescriptorBlock desc;

        // LogP block required, it is used in the Meylan LogP Model
        desc = new MeylanLogP();
        blocks.add(desc);

        // aggiungere blocco desc km??
        desc = new KmFactor();
        blocks.add(desc);


        return blocks;
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {

            // Calculate the whole Meylan LogP model
            LogPResult = LogPModel.Execute(CurMolecule, null, false);
            if (LogPResult.HasExperimental())
                LogPValueForModel = LogPResult.getExperimental();
            else
                LogPValueForModel = LogPResult.getMainResultValue();

            // Calculate the whole Meylan kM/half-life model
            kMResult = kMModel.Execute(CurMolecule, null, false);
            if (kMResult.HasExperimental())
                kMValueForModel = kMResult.getExperimental();
            else
                kMValueForModel = kMResult.getMainResultValue();

            Descriptors = new double[DescriptorsSize];
            Descriptors[0] = LogPResult.getMainResultValue();
            Descriptors[1] = kMResult.getMainResultValue();

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }



    @Override
    protected short CalculateModel() {

        if (LogPResult.getStatus() != InsilicoModelOutput.OUTPUT_OK)
            return MODEL_ERROR;

        if (kMResult.getStatus() != InsilicoModelOutput.OUTPUT_OK)
            return MODEL_ERROR;

        ArnotGobasBCFBAF arnotgobas = new ArnotGobasBCFBAF(LogPValueForModel, Math.pow(10, kMValueForModel));

        double MainResult = arnotgobas.getArnotLogBCFup();
        CurOutput.setMainResultValue(MainResult);

        String[] Res = new String[ResultsSize];

        Res[0] = String.valueOf(Format_2D.format(MainResult)); // BCF up
        double bufVal = Math.pow(10,MainResult);
        if (bufVal>1)
            Res[1] = String.valueOf(Math.round(bufVal)); // BCF not in log units (rounded to int)
        else
            Res[1] = String.valueOf(Format_2D.format(bufVal)); // BCF not in log units

        Res[2] = String.valueOf(Format_2D.format(arnotgobas.getArnotLogBCFlow())); // BCF low
        bufVal = Math.pow(10,arnotgobas.getArnotLogBCFlow());
        if (bufVal>1)
            Res[3] = String.valueOf(Math.round(bufVal)); // BCF not in log units (rounded to int)
        else
            Res[3] = String.valueOf(Format_2D.format(bufVal)); // BCF not in log units

        Res[4] = String.valueOf(Format_2D.format(arnotgobas.getArnotLogBCFmid())); // BCF mid
        bufVal = Math.pow(10,arnotgobas.getArnotLogBCFmid());
        if (bufVal>1)
            Res[5] = String.valueOf(Math.round(bufVal)); // BCF not in log units (rounded to int)
        else
            Res[5] = String.valueOf(Format_2D.format(bufVal)); // BCF not in log units

        Res[6] = String.valueOf(Format_2D.format(LogPValueForModel)); // LogP
        String Rel = "n.a.";
        if (LogPResult.HasExperimental())
            Rel = "Experimental";
        else
            switch (LogPResult.getADI().GetAssessmentClass()) {
                case ADIndex.INDEX_LOW:
                    Rel = "Low";
                    break;
                case ADIndex.INDEX_MEDIUM:
                    Rel = "Moderate";
                    break;
                case ADIndex.INDEX_HIGH:
                    Rel = "Good";
                    break;
            }
        Res[7] = Rel; // Reliability of logp prediction

        Res[8] = String.valueOf(Format_2D.format(kMValueForModel)); // LogP
        Rel = "n.a.";
        if (kMResult.HasExperimental())
            Rel = "Experimental";
        else
            switch (kMResult.getADI().GetAssessmentClass()) {
                case ADIndex.INDEX_LOW:
                    Rel = "Low";
                    break;
                case ADIndex.INDEX_MEDIUM:
                    Rel = "Moderate";
                    break;
                case ADIndex.INDEX_HIGH:
                    Rel = "Good";
                    break;
            }
        Res[9] = Rel; // Reliability of logp prediction

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.75);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.5);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }

        // LogP reliability
        ADIndexLogPReliability ADLogP = new ADIndexLogPReliability();
        ADLogP.SetReliability(LogPResult.getADI().GetAssessmentClass(), LogPResult.HasExperimental());
        CurOutput.addADIndex(ADLogP);

        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution;

        ADIndexADIMeylanBCF ADI = new ADIndexADIMeylanBCF();
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class),
                (iADIndex)ADLogP);
        CurOutput.setADI(ADI);

        // Reasoning items
        DescriptorAnalysisBCFMeylan DescLogp = new DescriptorAnalysisBCFMeylan(0);
        DescLogp.setDescriptorValue(Descriptors[0]);
        CurOutput.addReasoningItem(DescLogp);


        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

        double BCF_threshold_green = 2.7;
        double BCF_threshold_red = 3.3;

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0], " log(L/kg)");

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val < BCF_threshold_green)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val < BCF_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);

    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this, false, true);
        TSK.SerializeToFile(DatName);
    }
}
