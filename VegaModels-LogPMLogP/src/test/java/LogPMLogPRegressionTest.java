import insilico.core.model.InsilicoModel;
import insilico.logp_mlogp.ismLogPMLogP;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the LogPMLogP model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class LogPMLogPRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismLogPMLogP());
    }
}
