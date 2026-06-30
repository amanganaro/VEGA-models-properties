import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.meylanlogp.ismLogPMeylan;
import model.ModelExecutionTest;

public class LogPMeylanTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismLogPMeylan();
    }
}
