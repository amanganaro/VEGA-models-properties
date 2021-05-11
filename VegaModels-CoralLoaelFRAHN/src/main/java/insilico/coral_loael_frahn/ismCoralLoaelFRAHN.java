package insilico.coral_loael_frahn;

import com.hp.hpl.jena.assembler.Mode;
import insilico.coral.CoralModel;
import insilico.coral.models.loael.CoralLOAEL;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;

import java.util.ArrayList;

public class ismCoralLoaelFRAHN extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_coral_loael_frahn.xml";

    private CoralModel LoaelFrahn;

    public ismCoralLoaelFRAHN() throws InitFailureException{
        super(ModelData);

        try {
            LoaelFrahn = new CoralLOAEL();
        } catch (Exception ex){
            throw new InitFailureException("Unable to init CORAL model");
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];

        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];

        // check on result names and units
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted LOEL [log(mg/Kg body weight/day)]";
        this.ResultsName[1] = "Predicted LOEL [mg/Kg body weight/day]";
//        this.ResultsName[2] = "Molecules used for prediction";
//        this.ResultsName[3] = "Experimental LOEL [mg/Kg body weight/day]";


        // Define AD items
        this.ADItemsName = new String[0];


    }

    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {

        return new ArrayList<>();
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
        double Prediction;
        try {
            Prediction = LoaelFrahn.Predict(this.CurMolecule.GetSMILES());
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_3D.format(Prediction)); // -Log(mg/kg)
        double ConvertedValue = Math.pow(10, -1.0 * Prediction);
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // mg/kg
        else if (ConvertedValue>0.0001)
            Res[1] = Format_4D.format(ConvertedValue); // mg/kg
        else
            Res[1] = Format_6D.format(ConvertedValue); // mg/kg
        Res[2] = "-";
        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }

    @Override
    protected short CalculateAD() {
        return InsilicoModel.AD_ERROR;
    }

    @Override
    protected void CalculateAssessment() {
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // No assessment color for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
    }
}
