import insilico.core.model.InsilicoModel;
import insilico.henryslaw.ismHenrysLawOpera;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the HenrysLawOpera model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class HenrysLawOperaRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismHenrysLawOpera());
    }
}
