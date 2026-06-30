import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.watersolubility.ismWaterSolubilityIRFMN;
import model.ModelExecutionTest;

public class WaterSolubilityTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismWaterSolubilityIRFMN();
    }
}
