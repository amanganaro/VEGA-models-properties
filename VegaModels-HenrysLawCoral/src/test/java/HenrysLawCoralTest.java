import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.henryslawcoral.ismHenrysLawCoral;
import model.ModelExecutionTest;

public class HenrysLawCoralTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException {
        return new ismHenrysLawCoral();
    }
}
