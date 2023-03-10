package insilico.skin_irritation;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.ADCheckSA;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.skin_irritation.descriptors.EmbeddedDescriptors;
import insilico.skin_irritation.runner.nrNetwork;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.DataInputStream;
import java.net.URL;

public class ismSkinIrritation extends InsilicoModel {

    private static final Logger log = LogManager.getLogger(ismSkinIrritation.class);

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_irritation.xml";
    private static final String NNData = "/data/skinIrritation.nn";

    public ismSkinIrritation() throws InitFailureException {
        super(ModelData);
        this.DescriptorsSize = 44;
        this.DescriptorsNames = new String[DescriptorsSize];

        DescriptorsNames[0] = "SM1_Dz(i)";
        DescriptorsNames[1] = "Eig09_EA(ri)";
        DescriptorsNames[2] = "NdsCH";
        DescriptorsNames[3] = "X3Av";
        DescriptorsNames[4] = "NssS";
        DescriptorsNames[5] = "nRNR2";
        DescriptorsNames[6] = "CATS2D_01_DA";
        DescriptorsNames[7] = "IC3";
        DescriptorsNames[8] = "BB_SA44";
        DescriptorsNames[9] = "nArCHO";
        DescriptorsNames[10] = "nRCO";
        DescriptorsNames[11] = "SM03_AEA(dm)";
        DescriptorsNames[12] = "B03[N-O]";
        DescriptorsNames[13] = "B04[C-N]";
        DescriptorsNames[14] = "Eig05_EA(dm)";
        DescriptorsNames[15] = "C-026";
        DescriptorsNames[16] = "nCrs";
        DescriptorsNames[17] = "SpMaxA_D/Dt";
        DescriptorsNames[18] = "H-052";
        DescriptorsNames[19] = "nN";
        DescriptorsNames[20] = "CATS2D_09_AL";
        DescriptorsNames[21] = "Eig12_EA(dm)";
        DescriptorsNames[22] = "Eig15_EA(dm)";
        DescriptorsNames[23] = "CATS2D_08_NL";
        DescriptorsNames[24] = "SpMin6_Bh(s)";
        DescriptorsNames[25] = "SpMin6_Bh(i)";
        DescriptorsNames[26] = "piPC10";
        DescriptorsNames[27] = "B06[O-O]";
        DescriptorsNames[28] = "Eig09_AEA(ri)";
        DescriptorsNames[29] = "SpMAD_AEA(dm)";
        DescriptorsNames[30] = "B05[O-Cl]";
        DescriptorsNames[31] = "nR5";
        DescriptorsNames[32] = "nArOH";
        DescriptorsNames[33] = "B04[Cl-Cl]";
        DescriptorsNames[34] = "IC1";
        DescriptorsNames[35] = "MPC09";
        DescriptorsNames[36] = "B09[O-O]";
        DescriptorsNames[37] = "CATS2D_08_PL";
        DescriptorsNames[38] = "F08[C-N]";
        DescriptorsNames[39] = "GATS5s";
        DescriptorsNames[40] = "F03[C-Cl]";
        DescriptorsNames[41] = "H-051";
        DescriptorsNames[42] = "B06[C-N]";
        DescriptorsNames[43] = "GATS1e";
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Skin Irritation class";
    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors EmbeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            Descriptors[0] = EmbeddedDescriptors.SM1_Dzi;
            Descriptors[1] = EmbeddedDescriptors.Eig09_EAri;
            Descriptors[2] = EmbeddedDescriptors.NdsCH;
            Descriptors[3] = EmbeddedDescriptors.X3Av;
            Descriptors[4] = EmbeddedDescriptors.NssS;
            Descriptors[5] = EmbeddedDescriptors.nRNR2;
            Descriptors[6] = EmbeddedDescriptors.CATS2D_01_DA;
            Descriptors[7] = EmbeddedDescriptors.IC3;
            Descriptors[8] = EmbeddedDescriptors.BB_SA44;
            Descriptors[9] = EmbeddedDescriptors.nArCHO;
            Descriptors[10] = EmbeddedDescriptors.nRCO;
            Descriptors[11] = EmbeddedDescriptors.SM03_AEAdm;
            Descriptors[12] = EmbeddedDescriptors.B03_N_O;
            Descriptors[13] = EmbeddedDescriptors.B04_C_N;
            Descriptors[14] = EmbeddedDescriptors.Eig05_EAdm;
            Descriptors[15] = EmbeddedDescriptors.C_026;
            Descriptors[16] = EmbeddedDescriptors.nCrs;
            Descriptors[17] = EmbeddedDescriptors.SpMaxA_D_Dt;
            Descriptors[18] = EmbeddedDescriptors.H_052;
            Descriptors[19] = EmbeddedDescriptors.nN;
            Descriptors[20] = EmbeddedDescriptors.CATS2D_09_AL;
            Descriptors[21] = EmbeddedDescriptors.Eig12_EAdm;
            Descriptors[22] = EmbeddedDescriptors.Eig15_EAdm;
            Descriptors[23] = EmbeddedDescriptors.CATS2D_08_NL;
            Descriptors[24] = EmbeddedDescriptors.SpMin6_Bhs;
            Descriptors[25] = EmbeddedDescriptors.SpMin6_Bhi;
            Descriptors[26] = EmbeddedDescriptors.piPC10;
            Descriptors[27] = EmbeddedDescriptors.B06_O_O;
            Descriptors[28] = EmbeddedDescriptors.Eig09_AEAri;
            Descriptors[29] = EmbeddedDescriptors.SpMAD_AEAdm;
            Descriptors[30] = EmbeddedDescriptors.B05_O_Cl;
            Descriptors[31] = EmbeddedDescriptors.nR5;
            Descriptors[32] = EmbeddedDescriptors.nArOH;
            Descriptors[33] = EmbeddedDescriptors.B04_Cl_Cl;
            Descriptors[34] = EmbeddedDescriptors.IC1;
            Descriptors[35] = EmbeddedDescriptors.MPC09;
            Descriptors[36] = EmbeddedDescriptors.B09_O_O;
            Descriptors[37] = EmbeddedDescriptors.CATS2D_08_PL;
            Descriptors[38] = EmbeddedDescriptors.F08_C_N;
            Descriptors[39] = EmbeddedDescriptors.GATS5s;
            Descriptors[40] = EmbeddedDescriptors.F03_C_Cl;
            Descriptors[41] = EmbeddedDescriptors.H_051;
            Descriptors[42] = EmbeddedDescriptors.B06_C_N;
            Descriptors[43] = EmbeddedDescriptors.GATS1e;

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

            Prediction = nn.Calculate(Descriptors, true);
            CurOutput.setMainResultValue(Prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = this.GetTrainingSet().getClassLabel(Prediction);
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;

        } catch (Exception ex) {
            return MODEL_ERROR;
        }
    }

    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.AddMappingToNegativeValue(-1);
        adq.setMoleculesForIndexSize(3);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.6);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.9, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.9, 0.5);
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
        ADI.SetThresholds(0.9, 0.65);
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
