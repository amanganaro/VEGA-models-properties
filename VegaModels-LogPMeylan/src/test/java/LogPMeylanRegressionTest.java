import insilico.core.model.InsilicoModel;
import insilico.meylanlogp.ismLogPMeylan;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the LogPMeylan model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class LogPMeylanRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismLogPMeylan());
    }
}
