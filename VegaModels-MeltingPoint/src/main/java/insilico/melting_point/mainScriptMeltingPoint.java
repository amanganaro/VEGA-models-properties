package insilico.melting_point;

import insilico.core.exception.InitFailureException;
import utils.ModelsDeployment;

public class mainScriptMeltingPoint {

    public static void main (String[] args) throws InitFailureException {

        ismMeltingPoint model = new ismMeltingPoint();
        ModelsDeployment.BuildDataset(model, "out_ts");
    }
}
