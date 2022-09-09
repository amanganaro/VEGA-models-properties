import insilico.aromatase_irfmn.ismAromataseIRFMN;
import insilico.core.model.InsilicoModel;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class mainScriptAromatase {
    private static final Logger log = LogManager.getLogger(mainScriptAromatase.class);

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismAromataseIRFMN();
//        ModelsDeployment.BuildDataset(model, "out_ts");
//        File sourceFile = new File("out_ts/" + model.getInfo().getTrainingSetURL() + "/" + model.getInfo().getTrainingSetURL().split("/data/")[1]);
//        File destinationFile = new File("VegaModels-AromataseIRFMN\\src\\main\\resources\\data\\ts_aromatase_irfmn.dat");
//        try {
//            Files.move(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
//        } catch (Exception ex) {
//            log.warn(ex.getMessage());
//        }

//        List<String> smilesList = new ArrayList<>();
////        smilesList.add("O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)C");
////        smilesList.add("O=S(=O)(N)c1cc2c(cc1C(F)(F)F)NC(NS2(=O)(=O))Cc3ccccc3");
////        smilesList.add("O=[N+]([O-])c1cc(c(O)c(c1)C(C)(C)C)[N+](=O)[O-]");
////        smilesList.add("O=S(=O)(C(C)(C)C)C(C)(C)C");
//
//        smilesList.add("ClC1=C(C=CC([C@@H](O)[C@H]2NC([C@H](NC)CC(C)C)=O)=C1)OC3=C(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O[C@H]5C[C@@](NCC6=CC=C(C(NCCC[N+](C)(C)CC(NCCCCCCCCCC)=O)=O)C=C6)(C)[C@H](O)[C@H](C)O5)C(OC7=CC=C([C@@H](O)[C@H](NC8=O)C(N[C@@H](C9=CC(O)=CC(O)=C9%10)C(O)=O)=O)C=C7Cl)=CC([C@H](C(N[C@@H]8C%11=CC=C(O)C%10=C%11)=O)NC([C@H](CC(N)=O)NC2=O)=O)=C3");
//        smilesList.add("O=C1C23C(=NN1c4ccc(S(=O)(=O)O[Co-5]56728(C9(C(=O)N(c%10ccc(S(=O)(=O)O5)cc%10)N=C9C)NN=C%11C(=O)C=CC(S(=O)(=O)Nc%12ccc(N=NC7%13C(=O)N(c%14ccc(S(=O)(=O)O6)cc%14)N=C%13C)cc%12)=C%11)C%15(C(=O)N(c%16ccc(S(=O)(=O)O8)cc%16)N=C%15C)N=Nc%17ccc(cc%17)NS(=O)(=O)C%18=CC(C(=O)C=C%18)=NN3)cc4)C");
//        smilesList.add("C(=O)(N1[C@H](C(=O)O)CCC1)[C@@H](NC(=O)[C@H](C(C)C)NC(=O)[C@@H](NC(=O)[C@H](C(C)C)NC(=O)[C@H]2N(C(=O)[C@H](NC(=O)[C@H](NC(=O)[C@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@H](C(C)C)NC(=O)[C@H]3N(C(=O)[C@@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(CO)N)Cc4ccc(O)cc4)CO)CCSC)CCC(=O)O)Cc5cnc[nH]5)Cc6ccccc6)CCCN=C(N)N)Cc7c8c(cccc8)[nH]c7)CCCCN)CCC3)CCCCN)CCCCN)CCCN=C(N)N)CCCN=C(N)N)CCC2)CCCCN)Cc9ccc(O)cc9");
//
//
//        for(String smiles : smilesList) {
//            InsilicoModelOutput out = model.Execute(SmilesMolecule.Convert(smiles));
//            for(int i = 0; i < model.GetResultsName().length; i++){
//                System.out.println(smiles);
//                System.out.println(model.GetResultsName()[i] + " | " + out.getResults()[i]);
//            }
//        }

    }
}
