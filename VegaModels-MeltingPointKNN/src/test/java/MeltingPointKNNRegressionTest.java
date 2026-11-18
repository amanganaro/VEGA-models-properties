import insilico.core.model.InsilicoModel;
import insilico.melting_point_knn.ismMeltingPointKnn;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the MeltingPointKNN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class MeltingPointKNNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismMeltingPointKnn());
    }
}
