import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.watersolubility.ismWaterSolubilityIRFMN;

public class main {

    public static void main(String[]  args) throws Exception {

        InsilicoModel model = new ismWaterSolubilityIRFMN();
//        ModelsDeployment.BuildDataset(model, "out_ts");
        InsilicoModelOutput output = model.Execute(SmilesMolecule.Convert("C1CC1"));
        for(int i = 0; i < model.GetResultsSize(); i++) {
            System.out.println(model.GetResultsName()[i] + " | " + output.getResults()[i]);
        }
//        List<String> smilesList = new ArrayList<>();
//        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
//        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
//        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
//        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");



//        DescriptorBlock descriptorBlock = new MeylanLogP();
//        EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();

//        for(String smiles : smilesList) {
//            InsilicoMolecule mol = SmilesMolecule.Convert(smiles);
//            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors();
//            embeddedDescriptors.CalculateAllEmbeddedDescriptors(mol);
//            descriptorBlock.Calculate(mol);
//            System.out.println("LogP for " + smiles + ": " + embeddedDescriptors.BEH2p);
//        }

    }
}
