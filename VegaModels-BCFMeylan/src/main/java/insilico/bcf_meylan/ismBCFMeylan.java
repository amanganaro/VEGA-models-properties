package insilico.bcf_meylan;

import insilico.bcf_meylan.ad.ADIndexADIMeylanBCF;
import insilico.bcf_meylan.ad.ADIndexLogPReliability;
import insilico.bcf_meylan.ad.ADIndexMeylanBCFRange;
import insilico.bcf_meylan.descriptors.DescriptorAnalysisBCFMeylan;
import insilico.bcf_meylan.descriptors.DescriptorMW;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import insilico.meylanlogp.ismLogPMeylan;

/**
 *
 * @author User
 */
public class ismBCFMeylan extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_bcf_meylan.xml";
    
    private final BCFMeylanCorrectionsChecker BCFCorrection;
    private final BCFMeylanIonicChecker BCFIonic;
    private final ismLogPMeylan LogPModel;
    
    private InsilicoModelOutput LogPResult;
    private double LogPValueForModel;
    
    
    public ismBCFMeylan() 
            throws InitFailureException {
        super(ModelData);
        
        // Builds model objects
        LogPModel = new ismLogPMeylan();
        try {
            BCFCorrection = new BCFMeylanCorrectionsChecker();
            BCFIonic = new BCFMeylanIonicChecker();
        } catch (Exception ex) {
            throw new InitFailureException("unable to initialize internal SA checkers");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 2;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "LogP";
        this.DescriptorsNames[1] = "MW";
        
        // Defines results
        this.ResultsSize = 6;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted BCF [log(L/kg)]";
        this.ResultsName[1] = "Predicted BCF [L/kg]";
        this.ResultsName[2] = "Predicted LogP (Meylan/Kowwin)";
        this.ResultsName[3] = "Predicted LogP reliability";
        this.ResultsName[4] = "MW";
        this.ResultsName[5] = "Ionic compound";
        
        // Define AD items
        this.ADItemsName = new String[7];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexLogPReliability().GetIndexName();
        this.ADItemsName[5] = new ADIndexMeylanBCFRange().GetIndexName();
        this.ADItemsName[6] = new ADIndexACF().GetIndexName();
        
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
            
            // Calculate correction and ionic fragments
            BCFCorrection.SetLogP(LogPValueForModel);
            BCFCorrection.Calculate(CurMolecule);
            BCFIonic.Calculate(CurMolecule);
            
            Descriptors = new double[DescriptorsSize];
            DescriptorMW descriptorMW = new DescriptorMW(CurMolecule);
            Descriptors[0] = LogPResult.getMainResultValue();
            Descriptors[1] = descriptorMW.Mw;
            
            // MW in constitutional is given as a SCALED 
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            Descriptors[1] = Descriptors[1] * CarbonWeight;
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {

        if (LogPResult.getStatus() != InsilicoModelOutput.OUTPUT_OK)
            return MODEL_ERROR;
        
        double MainResult;

        if (BCFIonic.IsIonic()) {

            // IONIC compounds

            if (LogPValueForModel < 5.0)
                MainResult = 0.50;
            else if (LogPValueForModel < 6.0)
                MainResult = 0.75;
            else if (LogPValueForModel < 7.0)
                MainResult = 1.75;
            else if (LogPValueForModel < 9.0)
                MainResult = 1.00;
            else 
                MainResult = 0.50; // logP > 9.0

        } else {

            // NON-IONIC compounds

            if (LogPValueForModel > 7.0) {

                // Equation for logP > 7.0

                double Correction = BCFCorrection.GetCorrection();
                MainResult = (-0.49) * LogPValueForModel + 7.554 + Correction;

                // Minimum Tin and Mercury logBCF rule
                if (BCFCorrection.GetTinAndMercury())
                    if (MainResult < 2)
                        MainResult = 2;

                // Minimum logBCF for LogP>7
                if (MainResult < 0.5)
                    MainResult = 0.5;

            } else if (LogPValueForModel >= 1.0) {

                // Equation for logP [1.0 , 7.0]

                double Correction = BCFCorrection.GetCorrection();
                MainResult = (0.6598) * LogPValueForModel - 0.333 + Correction;

                // Minimum Tin and Mercury logBCF rule
                if (BCFCorrection.GetTinAndMercury())
                    if (MainResult < 2)
                        MainResult = 2;

            } else {

                // logP < 1.0

                MainResult = 0.50;  // fixed value
            }
        }        
        
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(MainResult)); // BCF
        double bufVal = Math.pow(10,MainResult);
        if (bufVal>1)
            Res[1] = String.valueOf(Math.round(bufVal)).toString(); // BCF not in log units (rounded to int)
        else
            Res[1] = String.valueOf(Format_2D.format(bufVal)); // BCF not in log units
        Res[2] = String.valueOf(Format_2D.format(LogPValueForModel)); // LogP
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
        Res[3] = Rel; // Reliability of logp prediction
        Res[4] = String.valueOf(Format_2D.format(Descriptors[1])); // MW
        Res[5] = BCFIonic.IsIonic() ? "yes" : "no";
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
        
        // Meylan BCF specific range check
        ADIndexMeylanBCFRange ADRange = new ADIndexMeylanBCFRange();
        ADRange.setDescriptors(Descriptors[0], Descriptors[1], BCFIonic.IsIonic());
        CurOutput.addADIndex(ADRange);
        
        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcContribution = ADRange.GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;

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
