package insilico.bcf_caesar;

import insilico.bcf_caesar.descriptors.BCFModel;
import insilico.bcf_caesar.descriptors.DescriptorAnalysisBCF;
import insilico.bcf_caesar.descriptors.EmbeddedDescriptors;
import insilico.bcf_caesar.descriptors.UncertaintyBCF;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SABCFCaesar;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.model.trainingset.TrainingSet;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class ismBCFCaesar extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_bcf_caesar.xml";
    
    private final BCFModel BCF;
    private final String BCFFullAlertSet;
    
    public ismBCFCaesar() 
            throws InitFailureException {
        super(ModelData);
        
        // Build model objects
        BCF = new BCFModel();
        
        // Build SA list
        BCFFullAlertSet = AlertEncoding.MergeAlertIds(new SABCFCaesar().getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 8;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "X0sol";
        this.DescriptorsNames[1] = "MATS5v";
        this.DescriptorsNames[2] = "GATS5v";
        this.DescriptorsNames[3] = "BEH2p";
        this.DescriptorsNames[4] = "AEige";
        this.DescriptorsNames[5] = "Cl-089";
        this.DescriptorsNames[6] = "MLogP";
        this.DescriptorsNames[7] = "ssCl";
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted BCF [log(L/kg)]";
        this.ResultsName[1] = "Predicted BCF [L/kg]";
        this.ResultsName[2] = "Predicted BCF from sub-model 1 (HM) [log(L/kg)]";
        this.ResultsName[3] = "Predicted BCF from sub-model 2 (GA) [log(L/kg)]";
        this.ResultsName[4] = "Predicted LogP (MLogP)";
        
        // Define AD items
        this.ADItemsName = new String[6];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexRange().GetIndexName();
        this.ADItemsName[5] = new ADIndexACF().GetIndexName();
        
    }
    

    
    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
        blocks.add(InsilicoConstants.SA_BLOCK_BCF_CAESAR);
        return blocks;
    }

    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList())
            if (AlertEncoding.ContainsAlert(BCFFullAlertSet, a.getId())) {
                FoundSAs.add((Alert)a.clone());
            }
        return FoundSAs;
    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.X0sol;
            Descriptors[1] = embeddedDescriptors.MATS5v;
            Descriptors[2] = embeddedDescriptors.GATS5v;
            Descriptors[3] = embeddedDescriptors.BEH2p;
            Descriptors[4] = embeddedDescriptors.AEige;
            Descriptors[5] = embeddedDescriptors.Cl_089;
            Descriptors[6] = embeddedDescriptors.MLogP;
            Descriptors[7] = embeddedDescriptors.ssCl;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        if (!BCF.Calculate(Descriptors)) 
            return MODEL_ERROR;
        
        CurOutput.setMainResultValue(BCF.getPrediction());
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(BCF.getPrediction())); // BCF
        double bufVal = Math.pow(10,BCF.getPrediction());
        if (bufVal>1)
            Res[1] = String.valueOf(Math.round(bufVal)).toString(); // BCF not in log units (rounded to int)
        else
            Res[1] = String.valueOf(Format_2D.format(bufVal)); // BCF not in log units
        Res[2] = String.valueOf(Format_2D.format(BCF.getPredictionHM())); // BCF from HM model
        Res[3] = String.valueOf(Format_2D.format(BCF.getPredictionGA())); // BCF from GA model
        Res[4] = String.valueOf(Format_2D.format(Descriptors[6])); // MLogP
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
        
        // Sets Range check
        ADCheckDescriptorRange adrc = new ADCheckDescriptorRange();
        if (!adrc.Calculate(TS, Descriptors, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.75, 1.0, 0.85, 0.75);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);
        
        // Sets SA
        try {
            ADCheckSA ADSA = new ADCheckSA(TS);
            if (!ADSA.Calculate(CurMolecule, GetCalculatedAlert(), CurOutput, adq.GetSimilarMolecules()))
                return InsilicoModel.AD_ERROR;
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
        // Reasoning items
        DescriptorAnalysisBCF DescLogp = new DescriptorAnalysisBCF(6);
        DescLogp.setDescriptorValue(Descriptors[6]);
        CurOutput.addReasoningItem(DescLogp);
        
        UncertaintyBCF UncBCF = new UncertaintyBCF();
        UncBCF.Evaluate(CurOutput.getMainResultValue(), ADI.GetIndexValue());
        CurOutput.addReasoningItem(UncBCF);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        double BCF_threshold_green = 2.7;
        double BCF_threshold_red = 3.3;
    
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0], "log(L/kg)");

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val < BCF_threshold_green)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val < BCF_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);

        // Additional assessment when some SAs are found
        try {
            AlertList curSA = GetCalculatedAlert();
            if (curSA.size() > 0) {
                for (Alert SA : curSA.getSAList() )
                    if (SA.getBoolProperty(InsilicoConstants.KEY_ALERT_BCF_CAESAR_OUTLIER)) {
                        CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() +
                                "\nWarning: the prediction may be not fully reliable due to the presence of one or more fragments related to model outliers.");
                        break;
                    }
                String SAs = "\nThe following relevant fragments have been found: " +
                        ModelUtilities.BuildSANameList(curSA.getSAList());
                CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() + SAs);
            }
        } catch (CloneNotSupportedException e) {
            log.warn("unable to add SAs in the assessment");
        }
        
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
