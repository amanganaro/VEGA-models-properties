package insilico.tpo_oberon;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.trainingset.TrainingSet;
import insilico.tpo_oberon.descriptors.EmbeddedDescriptors;
import insilico.tpo_oberon.knn.KnnAlgorithm;

public class ismTpoOberon extends InsilicoModel {

    private static final long serialVersionUID = 1L;
    private static final String ModelData = "/data/model_tpo_oberon.xml";

    private  KnnAlgorithm knn;

    public ismTpoOberon() throws InitFailureException {
        super(ModelData);

        knn = null;

        this.DescriptorsSize = 20;
        DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "GATS1e";
        this.DescriptorsNames[1] = "nArOH";
        this.DescriptorsNames[2] = "CATS2D_02_DL";
        this.DescriptorsNames[3] = "MATS1e";
        this.DescriptorsNames[4] = "MATS1s";
        this.DescriptorsNames[5] = "C-026";
        this.DescriptorsNames[6] = "CATS2D_03_DL";
        this.DescriptorsNames[7] = "B10[C-C]";
        this.DescriptorsNames[8] = "MATS1p";
        this.DescriptorsNames[9] = "nCb-";
        this.DescriptorsNames[10] = "nX";
        this.DescriptorsNames[11] = "Uc";
        this.DescriptorsNames[12] = "P_VSA_i_1";
        this.DescriptorsNames[13] = "SpMAD_B(v)";
        this.DescriptorsNames[14] = "nCbH";
        this.DescriptorsNames[15] = "GATS1s";
        this.DescriptorsNames[16] = "MATS1m";
        this.DescriptorsNames[17] = "MLOGP";
        this.DescriptorsNames[18] = "SpMax2_Bh(s)";
        this.DescriptorsNames[19] = "Eta_C_A";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted TPO";

    }

    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {


        try {
            Descriptors = new double[DescriptorsSize];
            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            this.Descriptors[0] = embeddedDescriptors.GATS1e;
            this.Descriptors[1] = embeddedDescriptors.nArOH;
            this.Descriptors[2] = embeddedDescriptors.CATS2D_02_DL;
            this.Descriptors[3] = embeddedDescriptors.MATS1e;
            this.Descriptors[4] = embeddedDescriptors.MATS1s;
            this.Descriptors[5] = embeddedDescriptors.C_026;
            this.Descriptors[6] = embeddedDescriptors.CATS2D_03_DL;
            this.Descriptors[7] = embeddedDescriptors.B10_C_C;
            this.Descriptors[8] = embeddedDescriptors.MATS1p;
            this.Descriptors[9] = embeddedDescriptors.nCb_;
            this.Descriptors[10] = embeddedDescriptors.nX;
            this.Descriptors[11] = embeddedDescriptors.Uc;
            this.Descriptors[12] = embeddedDescriptors.P_VSA_i_1;
            this.Descriptors[13] = embeddedDescriptors.SpMAD_B_v_;
            this.Descriptors[14] = embeddedDescriptors.nCbH;
            this.Descriptors[15] = embeddedDescriptors.GATS1s;
            this.Descriptors[16] = embeddedDescriptors.MATS1m;
            this.Descriptors[17] = embeddedDescriptors.MLOGP;
            this.Descriptors[18] = embeddedDescriptors.SpMax2_Bh_s;
            this.Descriptors[19] = embeddedDescriptors.Eta_C_A;

        } catch (Throwable e){
            return DESCRIPTORS_ERROR;
        }
        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {

        if (knn == null) {
            try {
                knn = new KnnAlgorithm(this, this.KnnSkipExperimental);
                knn.setkNeighbours(5);
            } catch (Exception e) {
                e.printStackTrace();
                return MODEL_ERROR;
            }
        }

        try {
            int prediction = knn.calculatePrediction(CurMolecule, Descriptors);

            CurOutput.setMainResultValue(prediction);

            String[] Res = new String[ResultsSize];
            Res[0] = this.GetTrainingSet().getClassLabel(prediction);
            CurOutput.setResults(Res);

            return MODEL_CALCULATED;
        } catch (InvalidMoleculeException | GenericFailureException e) {
            e.printStackTrace();
            return MODEL_ERROR;
        }
    }

    @Override
    protected short CalculateAD() {
        return 0;
    }

    @Override
    protected void CalculateAssessment() {

    }

    @Override
    public void ProcessTrainingSet() throws Exception {
        this.setSkipADandTSLoading(false);
        TrainingSet TSK = new TrainingSet();
        String TSPath = this.getInfo().getTrainingSetURL();
        String[] buf = TSPath.split("/");
        String DatName = buf[buf.length-1];
        TSPath = TSPath.substring(0, TSPath.length()-3) + "txt";
        TSK.Build(TSPath, this);
        TSK.SerializeToFile(DatName);
    }
}
