package insilico.skin_sensitization_toxtree;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.slf4j.Slf4j;

@Slf4j
public class ismSkinSensitization extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_skin_sensitization.xml";

    private final static boolean CALCULATE_AD = false;

    private ToxTreeSkinClassification ToxTreeSkin;


    public ismSkinSensitization() throws InitFailureException {
        super(ModelData);

        ToxTreeSkin = new ToxTreeSkinClassification();

        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Tox Tree classification prediction";
        this.ResultsName[1] = "Tox Tree classification name";
        this.ResultsName[2] = "Tox Tree classification description";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();



    }
    @Override
    protected short CalculateDescriptors(DescriptorsEngine descriptorsEngine) {
        try {

            Descriptors = new double[DescriptorsSize];

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }

    @Override
    protected short CalculateModel() {
        ToxTreeSkinClassification.SkinClass oToxTree;

        try {
            oToxTree = ToxTreeSkin.Calculate(CurMolecule);
        } catch (GenericFailureException e) {
            return MODEL_ERROR;
        }


        CurOutput.setMainResultValue(oToxTree.getId());
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(oToxTree.getId());
        Res[1] = oToxTree.getName();
        Res[2] = oToxTree.getDescription();

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;


    }

    @Override
    protected short CalculateAD() {
        // for now - AD not calculated
        return 0;
    }

    @Override
    protected void CalculateAssessment() { }
}
