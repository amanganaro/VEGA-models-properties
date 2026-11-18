import insilico.core.model.InsilicoModel;
import insilico.vapour_pressure.ismVapourPressure;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the VapourPressure model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class VapourPressureRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismVapourPressure());
    }
}
