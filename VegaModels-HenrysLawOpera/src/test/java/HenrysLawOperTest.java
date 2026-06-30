import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.henryslaw.ismHenrysLawOpera;
import model.ModelExecutionTest;

public class HenrysLawOperTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismHenrysLawOpera();
    }
}
