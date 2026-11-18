import insilico.core.model.InsilicoModel;
import insilico.melting_point.ismMeltingPoint;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the MeltingPoint model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class MeltingPointRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismMeltingPoint());
    }
}
