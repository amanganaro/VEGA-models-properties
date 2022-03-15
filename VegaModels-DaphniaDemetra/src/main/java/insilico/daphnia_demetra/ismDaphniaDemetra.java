package insilico.daphnia_demetra;

import insilico.core.model.trainingset.TrainingSet;
import insilico.daphnia_demetra.descriptors.EmbeddedDescriptors;
import insilico.daphnia_demetra.descriptors.FragmentsCount;
import insilico.daphnia_demetra.descriptors.MLogPdrgw;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//@Slf4j
public class ismDaphniaDemetra extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_daphnia_demetra.xml";

    private final DaphniaPredictor Daphnia;


    // Constructor - Init Model
    public ismDaphniaDemetra() throws InitFailureException {

        super(ModelData);

        //Builds model objects
        Daphnia = new DaphniaPredictor();

        //Define descriptors size and names
        this.DescriptorsSize = 16;

        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "BEH1m";
        this.DescriptorsNames[1] = "Eig1p";
        this.DescriptorsNames[2] = "IC2";
        this.DescriptorsNames[3] = "IDE";
        this.DescriptorsNames[4] = "MLogP";
        this.DescriptorsNames[5] = "Mp";
        this.DescriptorsNames[6] = "MW";
        this.DescriptorsNames[7] = "nHAcc";
        this.DescriptorsNames[8] = "nArNR2";
        this.DescriptorsNames[9] = "nP";
        this.DescriptorsNames[10] = "O_057";
        this.DescriptorsNames[11] = "O_060";
        this.DescriptorsNames[12] = "S_107";
        this.DescriptorsNames[13] = "SRW5";
        this.DescriptorsNames[14] = "T(F..Cl)";
        this.DescriptorsNames[15] = "WA";

        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted toxicity [-log(mol/l)]";
        this.ResultsName[1] = "Predicted toxicity [mg/l]";
        this.ResultsName[2] = "Molecular Weight";
        this.ResultsName[3] = "Experimental value [mg/l]";

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
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine)  {

        try {

            // Custom descriptors
            FragmentsCount DaphniaFrags = new FragmentsCount();
            DaphniaFrags.Calculate(CurMolecule);

            MLogPdrgw DaphniaMLogP = new MLogPdrgw();
            DaphniaMLogP.Calculate(CurMolecule);


            // Embedded descriptors
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();
            embeddedDescriptors.CalculateAllEmbeddedDescriptors(CurMolecule);

            Descriptors = new double[DescriptorsSize];

            Descriptors[0] = embeddedDescriptors.BEH1m;
            Descriptors[1] = embeddedDescriptors.Eig1p;
            Descriptors[2] = embeddedDescriptors.IC2;
            Descriptors[3] = embeddedDescriptors.IDE;
            Descriptors[4] = DaphniaMLogP.LogP;
            Descriptors[5] = CurMolecule.GetBasicDescriptorByName("Mp").getValue();
//            Descriptors[5] = constitutionalBlock.GetByName("Mp").getValue();
//            Descriptors[6] = constitutionalBlock.GetByName("MW").getValue();
//            Descriptors[6] = CurMolecule.GetBasicDescriptorByName("MW").getValue();
            Descriptors[6] = embeddedDescriptors.Mw;

            Descriptors[7] = DaphniaFrags.nHAcc;
            Descriptors[8] = DaphniaFrags.nArNR2;
//            Descriptors[9] = constitutionalBlock.GetByName("nP").getValue();
            Descriptors[9] = CurMolecule.GetBasicDescriptorByName("nP").getValue();
            Descriptors[10] = DaphniaFrags.O_057;
            Descriptors[11] = DaphniaFrags.O_060;
            Descriptors[12] = DaphniaFrags.S_107;
            Descriptors[13] = embeddedDescriptors.SRW5;
            Descriptors[14] = embeddedDescriptors.T_F_CL;
            Descriptors[15] = embeddedDescriptors.WA;

            // MW (descriptor 6) in constitutional is given as a SCALED 
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            Descriptors[6] = Descriptors[6] * CarbonWeight;

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;

    }

    @Override
    protected short CalculateModel() {

        double Prediction = Daphnia.predict(Descriptors);

        // Conversion to -log(mol/L) - TO DO: check here!
        Prediction = -1 * Math.log10(Math.pow(10, (-1 * Prediction)) / 1000);

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // toxicity -log(mol/L)

        try {
            double ConvertedValue = Math.pow(10, (-1 * Prediction)) * 1000 * CurMolecule.GetMolecularWeight();
            if (ConvertedValue>1)
                Res[1] = Format_2D.format(ConvertedValue); // in mg/L
            else
                Res[1] = Format_4D.format(ConvertedValue); // in mg/L
            Res[2] = Format_2D.format(CurMolecule.GetMolecularWeight()); // MW
            Res[3] = "-";
        } catch (InvalidMoleculeException ex){
//            log.warn(ex.getMessage());
        }

 // Converted experimental - set after
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;


    }

    @Override
    protected short CalculateAD() {

        // AD indices calculation
        ADCheckIndicesQuantitative adQuantitative = new ADCheckIndicesQuantitative(TS);
        adQuantitative.setMoleculesForIndexSize(2);
        if(!adQuantitative.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Thresholds for AD indices
        try {
            ((ADIndexSimilarity) CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.75, 0.65);
            ((ADIndexAccuracy) CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.7);
            ((ADIndexConcordance) CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.7);
            ((ADIndexMaxError) CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.7);

        } catch (Throwable ex){
            return InsilicoModel.AD_ERROR;
        }

        // Sets range check
        ADCheckDescriptorRange adCheckDescriptorRange = new ADCheckDescriptorRange();
        if(!adCheckDescriptorRange.Calculate(TS, Descriptors, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets ACF check
        ADCheckACF adCheckACF = new ADCheckACF(TS);
        if (!adCheckACF.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcfContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
        double ADIValue = adQuantitative.getIndexADI() * acfContribution * rcfContribution;

        // todo ADI of model is the same as BCF CAESAR
        // NON E VERO - TO DO
        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.8, 0.65, 1.0, 0.8, 0.65);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class), CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);

        // Add transformed (mg/L) experimental if needed
        if (CurOutput.HasExperimental()){
            double ConvertedExp = Math.pow(10, (-1 * CurOutput.getExperimental())) * 1000 * Descriptors[6];
            if (ConvertedExp > 1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedExp);
            else CurOutput.getResults()[3] = Format_4D.format(ConvertedExp);
        }

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {

        // [mg/L]
        double LC_threshold_red = 1;
        double LC_threshold_orange = 10;
        double LC_threshold_yellow= -100;

        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set in [mg/L] if available
        String ADItemWarnings = ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());

        String Result = CurOutput.getResults()[1] + " mg/l";

        switch (CurOutput.getADI().GetAssessmentClass()) {
            case ADIndex.INDEX_LOW:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_LOW, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_LOW, Result, ADItemWarnings));
                break;
            case ADIndex.INDEX_MEDIUM:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_MEDIUM, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_MEDIUM, Result, ADItemWarnings));
                break;
            case ADIndex.INDEX_HIGH:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_HIGH, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_HIGH, Result));
                if (!ADItemWarnings.isEmpty())
                    CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() +
                            String.format(MessagesAD.ASSESS_LONG_ADD_ISSUES, ADItemWarnings));
                break;
        }

        if (CurOutput.HasExperimental()){
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L"));
        }

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        Val = Math.pow(10, (-1 * Val)) * 1000 * Descriptors[6]; // convert to mg/L
        if (Val < LC_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val < LC_threshold_orange)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        else if (Val < LC_threshold_yellow)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);

    }





}
