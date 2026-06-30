import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.vapour_pressure.ismVapourPressure;
import model.ModelExecutionTest;

public class VapourPressureTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismVapourPressure();
    }
}
