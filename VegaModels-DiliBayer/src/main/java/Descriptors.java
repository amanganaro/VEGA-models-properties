import insilico.core.exception.GenericFailureException;
import insilico.core.python.Communication;
import insilico.core.tools.utils.FileUtilities;
import insilico.core.tools.utils.GeneralUtilities;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Descriptors {

    Communication communication;

    private static final Logger log = LogManager.getLogger(Descriptors.class);

    public Descriptors(String smiles) throws GenericFailureException {
        String uh=System.getProperty("user.home");
        communication = new Communication();
        communication.setAdditionalEnvVariables(Map.of(
                "PATH", uh+"\\miniconda3\\Scripts\\;"+uh+"\\miniconda3\\;"+
                "C:\\Program Files\\Python313\\Scripts\\;C:\\Program Files\\Python313\\;"));
        FileUtilities.WriteByteArrayToFile("input.csv", smiles.getBytes());
    }

    public boolean calculateDescriptors() throws GenericFailureException, IOException, InterruptedException, URISyntaxException {

        boolean isEnvSet = configureCondaEnv();
        boolean result = false;
        if (isEnvSet) {
            result = communication.executeCommandInCondaEnv("cddd", "cddd",
                    "--input input.csv", "--output descriptors.csv", "--smiles_header smiles");
            File file = new File("input.csv");
            file.delete();
        }

        return result;
    }

    public String getCondaEnv(){
        return "cddd";
    }



    public boolean configureCondaEnv() throws InterruptedException, IOException, URISyntaxException {
        boolean isSet=communication.checkCondaEnv(getCondaEnv());
        if(!isSet){
            URL resource = Descriptors.class.getResource("python" + File.separator + getCondaEnv()+".yml");
            if(resource != null){
                Path pathToEnvFile= Paths.get(resource.toURI()).toAbsolutePath();
                isSet = communication.configureCondaEnv(getCondaEnv(), pathToEnvFile);
            }
            else{
                log.error("Cannot find file {}", getCondaEnv()+".yml");
                return isSet;
            }

            // add default model folder to put the model data into the directory of cddd conda env
            String destination = System.getProperty("user.home");
            String source = System.getProperty("user.dir");
            if (System.getProperty("os.name").startsWith("Windows")) {
                destination += "\\AppData\\Local\\cddd\\cddd\\default_model";
                source += "\\python\\default_model";
            }
            else {
                destination += "/.local/share/cddd/default_model";
                source += "/python/default_model";
            }

            if(isSet){
                isSet = FileUtilities.copyExternalData(source, destination);
            }else{
                //TODO gestire con un eccezione?
            }
        }

        return isSet;
    }

}
