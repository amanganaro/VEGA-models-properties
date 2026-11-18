import insilico.core.model.InsilicoModel;
import insilico.logp_alogp.ismLogPALogP;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the LogPALogP model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class LogPALogPRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismLogPALogP());
    }
}
