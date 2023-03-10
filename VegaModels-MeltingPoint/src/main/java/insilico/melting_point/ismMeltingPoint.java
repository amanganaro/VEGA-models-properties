package insilico.melting_point;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.melting_point.descriptors.EmbeddedDescriptors;
import insilico.melting_point.runner.nrNetwork;
import javassist.runtime.Desc;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.DataInputStream;
import java.net.URL;

public class ismMeltingPoint extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismMeltingPoint.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_melting_point.xml";
    private static final String NNData = "/data/meltingPoint.nn";

    public ismMeltingPoint() throws InitFailureException {
        super(ModelData);

        // Define no. of descriptors
        this.DescriptorsSize = 20;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "piPC01";
        DescriptorsNames[1] = "CATS2D_00_DD";
        DescriptorsNames[2] = "ATS2i";
        DescriptorsNames[3] = "P_VSA_LogP_3";
        DescriptorsNames[4] = "nCIC";
        DescriptorsNames[5] = "CIC3";
        DescriptorsNames[6] = "GATS1e";
        DescriptorsNames[7] = "BIC1";
        DescriptorsNames[8] = "ATSC1i";
        DescriptorsNames[9] = "nArCOOH";
        DescriptorsNames[10] = "NPerc";
        DescriptorsNames[11] = "B02[C-O]";
        DescriptorsNames[12] = "ATS3m";
        DescriptorsNames[13] = "nCconj";
        DescriptorsNames[14] = "F04[C-O]";
        DescriptorsNames[15] = "nROH";
        DescriptorsNames[16] = "CATS2D_04_PL";
        DescriptorsNames[17] = "S-106";
        DescriptorsNames[18] = "T(N..N)";
        DescriptorsNames[19] = "N-072";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Melting Point [°C]";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = embeddedDescriptors.piPC01;
            Descriptors[1] = embeddedDescriptors.CATS2D_00_DD;
            Descriptors[2] = embeddedDescriptors.ATS2i;
            Descriptors[3] = embeddedDescriptors.P_VSA_LogP_3;
            Descriptors[4] = embeddedDescriptors.nCIC;
            Descriptors[5] = embeddedDescriptors.CIC3;
            Descriptors[6] = embeddedDescriptors.GATS1e;
            Descriptors[7] = embeddedDescriptors.BIC1;
            Descriptors[8] = embeddedDescriptors.ATSC1i;
            Descriptors[9] = embeddedDescriptors.nArCOOH;
            Descriptors[10] = embeddedDescriptors.NPerc;
            Descriptors[11] = embeddedDescriptors.B02_C_O;
            Descriptors[12] = embeddedDescriptors.ATS3m;
            Descriptors[13] = embeddedDescriptors.nCconj;
            Descriptors[14] = embeddedDescriptors.F04_C_O;
            Descriptors[15] = embeddedDescriptors.nROH;
            Descriptors[16] = embeddedDescriptors.CATS2D_04_PL;
            Descriptors[17] = embeddedDescriptors.S_106;
            Descriptors[18] = embeddedDescriptors.T_N_N;
            Descriptors[19] = embeddedDescriptors.N_072;
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        double Prediction;
        try {
            DataInputStream in;
            URL nnURL = getClass().getResource(NNData);
            in = new DataInputStream(nnURL.openStream());
            nrNetwork nn = nrNetwork.ReadFromFile(in);

            Prediction = nn.Calculate(Descriptors,true);

            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = String.valueOf(Format_2D.format(Prediction));
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex){
            return MODEL_ERROR;
        }

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
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(50.0, 10.0);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(50.0, 10.0);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(50.0, 10.0);
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

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;

    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0], "°C");

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
