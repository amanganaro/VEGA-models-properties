import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.melting_point.ismMeltingPoint;
import model.ModelExecutionTest;

public class MeltingPointTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismMeltingPoint();
    }
}
