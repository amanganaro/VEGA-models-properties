package insilico.carcinogenicity_rat_male;

import insilico.carcinogenicity_rat_male.utils.ModelsDeployment;
import insilico.core.coral.CoralModel;
import insilico.core.coral.models.carcino.CoralMRCarcinogenicity;
import insilico.core.model.InsilicoModel;
import insilico.core.molecule.conversion.SmilesMolecule;

public class mainScript {

    public static void main(String[] args) throws Exception {
        InsilicoModel model = new ismCarcinogenicityRatMale();
//        ModelsDeployment.BuildDataset(model, "out_ts");

////        ModelsDeployment modelsDeployment = new ModelsDeployment().PrintDescriptor(model, "descriptors");
////        ModelsDeployment modelsDeployment = new ModelsDeployment().TestModelWithTrainingSet(model, "results_descriptors_embedded");
//
//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("c1(cc2c(cc1Cl)Oc1c(cc(c(c1)Cl)Cl)O2)Cl");
//        smilesList.add("C1(=O)CCN(N=O)C(=O)N1");
//        smilesList.add("C1=CCCN(C1)N=O");
//        smilesList.add("C1CCCN1NO");
//
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++)
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//        }
//        CoralModel coralModel = new CoralMRCarcinogenicity();
//        double val = coralModel.Predict("N(CCO)(N=O)CCCl");
//        double val_vega = coralModel.Predict(SmilesMolecule.Convert("N(CCO)(N=O)CCCl").GetSMILES());

//        return;
        ModelsDeployment.TestModelWithTrainingSet(model, "carcinogenicity_male");


    }

}
