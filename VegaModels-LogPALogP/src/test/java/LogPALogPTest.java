import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.logp_alogp.ismLogPALogP;
import model.ModelExecutionTest;

public class LogPALogPTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismLogPALogP();
    }
}
