package insilico.steroidogenesisedscreen.model;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.steroidogenesisedscreen.model.modelKNN.MolecularKNNClassifierFixed;
import insilico.steroidogenesisedscreen.model.modelRFpmml.RFResult;
import insilico.steroidogenesisedscreen.model.modelRFpmml.SteroidogenesisRFmodel;


public class SteroidogenesisConsModel {

    private final MolecularKNNClassifierFixed model_knn;
    private final SteroidogenesisRFmodel model_rf;

    public SteroidogenesisConsModel() throws InitFailureException {
        try {
            this.model_knn = new MolecularKNNClassifierFixed();
        } catch (Throwable e) {
            throw new InitFailureException("Unable to init KNN model - " + e.getMessage());
        }
        try {
            this.model_rf = new SteroidogenesisRFmodel();
        } catch (Throwable e) {
            throw new InitFailureException("Unable to init RF model - " + e.getMessage());
        }
    }


    public SteroidogenesisConsResult predict(InsilicoMolecule mol) {

        SteroidogenesisConsResult res = new SteroidogenesisConsResult();

        try {
            MolecularKNNClassifierFixed.Prediction pred_knn = this.model_knn.predict(mol.GetSMILES());
            if (pred_knn.prediction.equalsIgnoreCase("0"))
                res.KNN_result = 0;
            else if (pred_knn.prediction.equalsIgnoreCase("1"))
                res.KNN_result = 1;
            else
                res.KNN_result = Descriptor.MISSING_VALUE;
            res.KNN_similarity = pred_knn.similarity;
            res.KNN_n_positive = pred_knn.positives;
            res.KNN_n_negative = pred_knn.negatives;
        } catch (Throwable e) {
            // do nothing - knn values in results object already initialized to missing value
        }

        try {
            RFResult pred_rf = this.model_rf.predict(mol.getInputSMILES());
            res.RF_result = pred_rf.Prediction;
            res.RF_confidence = pred_rf.Confidence;
        } catch (Throwable e) {
            // do nothing - knn values in results object already initialized to missing value
        }

        // Combine results
        if  ( (res.RF_result == 1) && (res.KNN_result == 1) )
            res.Result = 1;
        else if  ( (res.RF_result == 0) && (res.KNN_result == 0) )
            res.Result = 0;
        else if  ( (res.RF_result == 1) && (res.KNN_result == 0) )
            res.Result = Descriptor.MISSING_VALUE;
        else if  ( (res.RF_result == 0) && (res.KNN_result == 1) )
            res.Result = Descriptor.MISSING_VALUE;
        else if  ( (res.RF_result == Descriptor.MISSING_VALUE) && (res.KNN_result != Descriptor.MISSING_VALUE) )
            res.Result = res.KNN_result;
        else if  ( (res.RF_result != Descriptor.MISSING_VALUE) && (res.KNN_result == Descriptor.MISSING_VALUE) )
            res.Result = res.RF_result;
        else
            res.Result = Descriptor.MISSING_VALUE;

        if (res.Result == Descriptor.MISSING_VALUE)
            res.Result = -1;

        return res;
    }

}
