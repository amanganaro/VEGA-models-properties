package insilico.ld50;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.ld50.utils.ModelsDeployment;

import java.io.*;
import java.net.MalformedURLException;

public class mainScriptLD50 {

    public static void main (String[] args) throws Exception {

        InsilicoModel model = new ismLD50();
        model.SetKnnSkipExperimental(false);
//
//        BufferedReader br = new BufferedReader(new FileReader("C:\\Users\\alber\\Desktop\\prova\\ld50.txt"));
//        String buf;
//        int idx = 0;
//        while ( (buf = br.readLine()) != null) {
//            idx++;
//            InsilicoMolecule m = SmilesMolecule.Convert(buf);
//            InsilicoModelOutput o = model.Execute(m);
//            String r = idx + "\t" + buf + "\t";
//            r += o.getResults()[0];
//            System.out.println(r);
//        }


        InsilicoMolecule m = SmilesMolecule.Convert("O=C(O)CC(c1ccc(cc1)Cl)CN");
        InsilicoModelOutput o = model.Execute(m);
        System.out.println(o.getResults()[0]);
        System.out.println(o.getResults()[1]);
        System.out.println(o.getResults()[2]);
        System.out.println(o.getAssessment());
        System.out.println(o.getAssessmentVerbose());
        if (1==1) return;

//        ModelsDeployment.BuildDataset(model, "out_ts");
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results", false);
        ModelsDeployment.TestModelWithTrainingSet(model, "ld50_knn_results_kode_dataset", true);

    }
}
