package insilico.tpo_oberon.knn;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.localization.StringSelectorCore;
import insilico.core.model.InsilicoModel;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.similarity.Similarity;
import insilico.core.similarity.SimilarityDescriptorsBuilder;
import insilico.core.tools.utils.GeneralUtilities;
import insilico.tpo_oberon.descriptors.EmbeddedDescriptors;
import insilico.tpo_oberon.descriptors.weights.DescriptorMLogP;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;


public class KnnAlgorithm {
    private static final Logger log = LogManager.getLogger(KnnAlgorithm.class);

    private final static int DESC_NUM = 20;

    private int kNeighbours;
    private static List<DatasetElement> datasetElementList;
    private final boolean skipExp;

    public KnnAlgorithm(InsilicoModel modelParent, boolean skipExp) throws Exception {
        this.skipExp = skipExp;
        this.kNeighbours = 5;
        initializeDescriptors(modelParent);
    }

    public void setkNeighbours(int kNeighbours) {
        this.kNeighbours = kNeighbours;
    }

    private void initializeDescriptors(InsilicoModel modelParent) throws Exception {

        ArrayList<String> SMI = new ArrayList<>();
        ArrayList<Integer> Exp = new ArrayList<>();
        ArrayList<Integer> Status = new ArrayList<>();

        iTrainingSet trainingSet = modelParent.GetTrainingSet();

        // check if TS is null (not available as model is in the deployment phase)
        if (trainingSet == null) {

            // read needed info directly from txt ts
            URL tsURL = this.getClass().getResource(modelParent.getInfo().getTrainingSetURL());
            DataInputStream in = new DataInputStream(tsURL.openStream());
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(in));
            String[] parsedString = bufferedReader.readLine().split("\t");
            if (parsedString.length < 5) {
                throw new GenericFailureException(StringSelectorCore.getString("trainingset_header_error"));
            }

            String string;
            while((string = bufferedReader.readLine()) != null) {
                parsedString = GeneralUtilities.TrimString(string).split("\t");
                SMI.add(parsedString[2]);
                if (parsedString[3].compareToIgnoreCase("Training") == 0) {
                    Status.add(1);
                } else if (parsedString[3].compareToIgnoreCase("Test") == 0) {
                    Status.add(2);
                } else {
                    Status.add(-1);
                }
                Exp.add(Integer.parseInt(parsedString[4]));
            }

        } else {

            for (int i=0; i<trainingSet.getMoleculesSize(); i++) {
                SMI.add(trainingSet.getSMILES(i));
                Status.add((int)trainingSet.getMoleculeSet(i));
                Exp.add((int)trainingSet.getExperimentalValue(i));
            }

        }

        datasetElementList = new ArrayList<>();
        SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = false;

        for (int dsIndex=0; dsIndex<SMI.size(); dsIndex++) {
            InsilicoMolecule curMolecule = SmilesMolecule.Convert(SMI.get(dsIndex));

            double[] descriptors = new double[DESC_NUM];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(curMolecule);
            descriptors[0] = embeddedDescriptors.GATS1e;
            descriptors[1] = embeddedDescriptors.nArOH;
            descriptors[2] = embeddedDescriptors.CATS2D_02_DL;
            descriptors[3] = embeddedDescriptors.MATS1e;
            descriptors[4] = embeddedDescriptors.MATS1s;
            descriptors[5] = embeddedDescriptors.C_026;
            descriptors[6] = embeddedDescriptors.CATS2D_03_DL;
            descriptors[7] = embeddedDescriptors.B10_C_C;
            descriptors[8] = embeddedDescriptors.MATS1p;
            descriptors[9] = embeddedDescriptors.nCb_;
            descriptors[10] = embeddedDescriptors.nX;
            descriptors[11] = embeddedDescriptors.Uc;
            descriptors[12] = embeddedDescriptors.P_VSA_i_1;
            descriptors[13] = embeddedDescriptors.SpMAD_B_v_;
            descriptors[14] = embeddedDescriptors.nCbH;
            descriptors[15] = embeddedDescriptors.GATS1s;
            descriptors[16] = embeddedDescriptors.MATS1m;
            descriptors[17] = embeddedDescriptors.MLOGP;
            descriptors[18] = embeddedDescriptors.SpMax2_Bh_s;
            descriptors[19] = embeddedDescriptors.Eta_C_A;

            datasetElementList.add(new DatasetElement(
                    curMolecule,
                    descriptors,
                    Exp.get(dsIndex),
                    Status.get(dsIndex)
            ));
            dsIndex++;
        }

        SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = true;



//        URL url = getClass().getResource("/data/dataset_descriptors.csv");
//
//        // calculate descriptors for all molecules in dataset
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){
//
//            br.readLine();
//            String line;
//
//            int dsIndex = 0;
//
//            while ((line = br.readLine()) !=null){
//
//                String[] curLine = line.split("\t");
//                String smiles = curLine[1];
//
//                SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = false;
//                InsilicoMolecule curMolecule = SmilesMolecule.Convert(smiles);
//
//                List<String> descriptorArrayList = Arrays.asList(curLine).subList(2, 2 + DESC_NUM);
//                double[] descriptors = new double[descriptorArrayList.size()];
//                for(int i = 0; i < descriptors.length; i++){
//                    descriptors[i] = Double.parseDouble(descriptorArrayList.get(i));
//                }
//
//                datasetElementList.add(new DatasetElement(
//                        curMolecule,
//                        descriptors,
//                        (int) trainingSet.getExperimentalValue(dsIndex),
//                        trainingSet.getMoleculeSet(dsIndex)
//                ));
//                dsIndex++;
//            }
//
//        } catch (GenericFailureException | IOException ex) {
//            System.out.println(ex.getMessage());
//        }
    }

    public Integer calculatePrediction(InsilicoMolecule curMolecule, double[] curMoleculeDescriptors) throws InvalidMoleculeException {

        List<DistanceClassification> distances = new ArrayList<>();

        try {
            // calculate distance for every dataset element
            for(DatasetElement curDatasetElement : datasetElementList) {

                double euclideanDistance = calculateDistance(curMoleculeDescriptors, curDatasetElement.getDescriptors());

                // if KNN descriptors-based distance is zero, check with VEGA similarity
                if (euclideanDistance == 0) {
                    SimilarityDescriptorsBuilder similarityDescriptorsBuilder = new SimilarityDescriptorsBuilder();
                    Similarity similarity = new Similarity();
                    double sim = similarity.CalculateExactMatches(
                            similarityDescriptorsBuilder.Calculate(curMolecule),
                            similarityDescriptorsBuilder.Calculate(curDatasetElement.getInputMolecule()),
                            curMolecule.GetStructure(),
                            curDatasetElement.getInputMolecule().GetStructure()
                    );

                    if (sim == 1.0) {
                        if (this.skipExp)
                            continue;
                        else return curDatasetElement.getExperimentalValue();
                    }
                }

                distances.add(new DistanceClassification(euclideanDistance, curDatasetElement.getExperimentalValue()));
            }

            distances.sort(Comparator.comparing(DistanceClassification::getDistance));

            List<DistanceClassification> selectedDistances = distances.subList(0, kNeighbours);

            double weight_0 = 0.0;
            double weight_1 = 0.0;

            for(DistanceClassification neigh : selectedDistances){
                if(neigh.getClassification() == 0)
                    weight_0 += neigh.getWeight();
                else
                    weight_1 += neigh.getWeight();
            }

            return weight_0 > weight_1 ? 0 : 1;

//            if(weight_0 > weight_1)
//                return 0;
//            else return 1;
        } catch (Exception ex){
            System.out.println(ex.getMessage());
            return -1;
        }
    }


    private double calculateDistance(double[] curMoleculeDescriptors, double[] curDatasetElementDescriptors){
        double distance = 0.0;

        for(int i = 0; i < curMoleculeDescriptors.length; i++){
            distance += Math.pow(curMoleculeDescriptors[i] - curDatasetElementDescriptors[i], 2);
        }
        return Math.sqrt(distance);
    }

}