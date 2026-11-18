import insilico.core.model.InsilicoModel;
import insilico.koc_opera.ismKocOpera;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the KocOpera model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class KocOperaRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismKocOpera());
    }
}
