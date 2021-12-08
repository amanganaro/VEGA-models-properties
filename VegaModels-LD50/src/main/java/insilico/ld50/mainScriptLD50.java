package insilico.ld50;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.ld50.utils.ModelsDeployment;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

public class mainScriptLD50 {

    public static void main (String[] args) throws InitFailureException, MalformedURLException, FileNotFoundException, GenericFailureException {

        InsilicoModel model = new ismLD50();
        InsilicoMolecule m = SmilesMolecule.Convert("C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)");
        InsilicoModelOutput o = model.Execute(m);
        System.out.println(o.getAssessment());
        System.out.println(o.getAssessmentVerbose());
        if (1==1) return;

//        ModelsDeployment.BuildDataset(model, "out_ts");
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results", false);
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results_kode_dataset", true);

    }
}
