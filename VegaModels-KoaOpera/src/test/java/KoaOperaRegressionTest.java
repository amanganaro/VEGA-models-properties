import insilico.core.model.InsilicoModel;
import insilico.koa_opera.ismKoaOpera;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the KoaOpera model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class KoaOperaRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismKoaOpera());
    }
}
