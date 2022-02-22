package insilico.tpo_oberon.knn;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.similarity.Similarity;
import insilico.core.similarity.SimilarityDescriptorsBuilder;
import lombok.extern.slf4j.Slf4j;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

@Slf4j
public class KnnAlgorithm {

    private int kNeighbours;
    private static List<DatasetElement> datasetElementList;
    private final double[] curMoleculeDescriptors;

    private final InsilicoMolecule curMolecule;
    private final boolean skipExp;

    public KnnAlgorithm(InsilicoMolecule insilicoMolecule, double[] curMoleculeDescriptors, iTrainingSet trainingSet, boolean skipExp) {
        this.curMolecule = insilicoMolecule;
        this.skipExp = skipExp;
        this.curMoleculeDescriptors = curMoleculeDescriptors;
        this.kNeighbours = 3;
        initializeDescriptors(trainingSet);
    }

    public void setkNeighbours(int kNeighbours) {
        this.kNeighbours = kNeighbours;
    }

    private void initializeDescriptors(iTrainingSet trainingSet) {

        datasetElementList = new ArrayList<>();
        URL url = getClass().getResource("/data/dataset_descriptors.csv");

        // calculate descriptors for all molecules in dataset
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(url.openStream())))){

            br.readLine();
            String line;

            int dsIndex = 0;

            while ((line = br.readLine()) !=null){

                String[] curLine = line.split("\t");
                String smiles = curLine[1];

                SmilesMolecule.EXCLUDE_DISCONNECTED_STRUCTURES = false;
                InsilicoMolecule curMolecule = SmilesMolecule.Convert(smiles);

                List<String> descriptorArrayList = Arrays.asList(curLine).subList(2, 2 + curMoleculeDescriptors.length);
                double[] descriptors = new double[descriptorArrayList.size()];
                for(int i = 0; i < descriptors.length; i++){
                    descriptors[i] = Double.parseDouble(descriptorArrayList.get(i));
                }

                datasetElementList.add(new DatasetElement(
                        curMolecule,
                        descriptors,
                        (int) trainingSet.getExperimentalValue(dsIndex),
                        trainingSet.getMoleculeSet(dsIndex)
                ));
                dsIndex++;
            }

        } catch (GenericFailureException | IOException ex) {
            System.out.println(ex.getMessage());
        }
    }

    public Integer calculatePrediction() throws InvalidMoleculeException {

        List<DistanceClassification> distances = new ArrayList<>();

        try {
            // calculate distance for every dataset element
            for(DatasetElement curDatasetElement : datasetElementList) {


                SimilarityDescriptorsBuilder similarityDescriptorsBuilder = new SimilarityDescriptorsBuilder();
                Similarity similarity = new Similarity();
                double sim = similarity.CalculateExactMatches(
                        similarityDescriptorsBuilder.Calculate(curMolecule),
                        similarityDescriptorsBuilder.Calculate(curDatasetElement.getInputMolecule()),
                        curMolecule.GetStructure(),
                        curDatasetElement.getInputMolecule().GetStructure()
                );

                if(sim == 1.0){
                    if (this.skipExp)
                        continue;
                    else return curDatasetElement.getExperimentalValue();
                }




                double euclideanDistance = calculateDistance(curMoleculeDescriptors, curDatasetElement.getDescriptors());
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