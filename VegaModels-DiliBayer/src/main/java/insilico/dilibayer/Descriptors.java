package insilico.dilibayer;

import insilico.core.exception.GenericFailureException;
import insilico.core.python.Communication;
import insilico.core.tools.utils.FileUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class Descriptors {

    Communication communication;
    private static final Logger log = LoggerFactory.getLogger(Descriptors.class);
    private final ismDiliBayer ismDiliBayer;
    private Path pathToExternalFolder;

    public Descriptors(ismDiliBayer model) {
        ismDiliBayer = model;

        String uh=System.getProperty("user.home");
        communication = new Communication();
        communication.setAdditionalEnvVariables(Map.of(
                "PATH", uh+"\\miniconda3\\Scripts\\;"+uh+"\\miniconda3\\;"+
                "C:\\Program Files\\Python313\\Scripts\\;C:\\Program Files\\Python313\\;"));

        if (System.getProperty("os.name").startsWith("Windows")) {
            pathToExternalFolder = Paths.get(System.getProperty("user.home"),
                    "\\AppData\\Local\\vega-models\\dili-bayer\\python\\").resolve("");
        }
        else {
            pathToExternalFolder = Paths.get(System.getProperty("user.home") ,"/.local/share/vega-models/dili-bayer/python/");
        }
    }

    public boolean calculateDescriptors() throws IOException, InterruptedException, URISyntaxException {
        boolean isEnvSet = ismDiliBayer.CHECK_SETUP ? configureCondaEnv() : true;
        boolean result = false;
        if (isEnvSet) {
            log.info("Start to calculate descriptors");
            result = communication.executeCommandInCondaEnv("cddd", "cddd",
                    "--input "+ismDiliBayer.getInputTempFile(), " --output "+ ismDiliBayer.getDescriptorsTempFile(),
                    " --smiles_header smiles");
            File f = new File(ismDiliBayer.getInputTempFile());
            f.delete();
            log.info("Finish to calculate descriptors");
        }
        return result;
    }

    public String getCondaEnv(){
        return "cddd";
    }


    /***
     * Move the .yml .whl and default_model folder into Local data folder to use that files to
     * setup the configuration. This is made to use external data instead the one in the project
     * @return
     * @throws InterruptedException
     * @throws IOException
     * @throws URISyntaxException
     */
    public boolean configureCondaEnv() throws InterruptedException, IOException, URISyntaxException {
        boolean isSet=communication.checkCondaEnv(getCondaEnv());
        if(!isSet){
            URL urlEnv = getClass().getResource("/python/"+getCondaEnv()+".yml");
            URL urlWheel = getClass().getResource("/python/cddd-1.2.3-py3-none-any.whl");
            URL urlModelDefaultFolder = getClass().getResource("/python/default_model");

            if(urlEnv!=null && urlWheel != null && urlModelDefaultFolder != null){
                FileUtilities.copyExternalData(Paths.get(urlEnv.toURI()).toString(),
                        pathToExternalFolder.toString());
                FileUtilities.copyExternalData(Paths.get(urlWheel.toURI()).toString(),
                        pathToExternalFolder.toAbsolutePath().toString());

                Path pathToEnvFile= Paths.get(urlEnv.toURI()).toAbsolutePath();
                isSet = communication.configureCondaEnv(getCondaEnv(), pathToEnvFile);
                if(isSet){
                    // add default model folder to put the model data into the directory of cddd conda env
                    String destination = System.getProperty("user.home");
                    if (System.getProperty("os.name").startsWith("Windows")) {
                        destination += "\\AppData\\Local\\cddd\\cddd\\default_model\\";
                    }
                    else {
                        destination += "/.local/share/cddd/default_model/";
                    }

                    isSet = FileUtilities.copyExternalData(Paths.get(urlModelDefaultFolder.toURI()).toString(),
                            destination);
                }else{
                    log.error("Error in set up conda environment {}",getCondaEnv());
                }
            }
            else{
                log.error("Some files to setup cddd conda environment {} are missing", getCondaEnv());
            }
        }

        log.info("Conda environment {} set up {}", getCondaEnv(), isSet ? "correctly": "failed");

        return isSet;
    }

}
