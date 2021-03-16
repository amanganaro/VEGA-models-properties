import insilico.core.model.InsilicoModel;
import insilico.core.molecule.acf.ACFItem;
import insilico.daphnia_demetra.ismDaphniaDemetra;
import utils.ModelsDeployment;

import java.util.ArrayList;

public class main {

    public static void main(String[] args) throws Exception {


        InsilicoModel model = new ismDaphniaDemetra();
//        ModelsDeployment.BuildDataset(models, "out_ts");
        for(ACFItem acf : model.GetTrainingSet().getACF().getList()){
            System.out.println(acf.getACF() + " : " + acf.getFrequency());
        }

    }
}
