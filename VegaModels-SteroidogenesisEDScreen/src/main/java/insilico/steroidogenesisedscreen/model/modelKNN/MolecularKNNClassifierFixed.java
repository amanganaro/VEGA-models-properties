package insilico.steroidogenesisedscreen.model.modelKNN;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.*;
import java.util.stream.Collectors;

public class MolecularKNNClassifierFixed {
    
    // Fixed parameters
    private static final int K = 5;
    private static final double MIN_SIMILARITY = 0.6;
    
    private static class MoleculeData {
        String ds, id, smiles, target;
        BitSet fingerprint;
        
        public MoleculeData(String ds, String id, String smiles, String target) {
            this.ds = ds;
            this.id = id;
            this.smiles = smiles;
            this.target = target;
            this.fingerprint = null; // at beginning all FPs are set to null
        }
    }
    
    public static class Prediction {
        public String ds, id, smiles, target, prediction;
        public double similarity;
        public int neighbors, positives, negatives;
        
        public Prediction(String ds, String id, String smiles, String target, String prediction, 
                         double similarity, int neighbors, int positives, int negatives) {
            this.ds = ds;
            this.id = id;
            this.smiles = smiles;
            this.target = target;
            this.prediction = prediction;
            this.similarity = similarity;
            this.neighbors = neighbors;
            this.positives = positives;
            this.negatives = negatives;
        }
    }
    
    private static class ClassificationStats {
        String dataset;
        int tp, fp, tn, fn, outOfDomain;
        double coverage, balancedAccuracy, accuracy, sensitivity, specificity, ppv, npv, mcc;
        
        public ClassificationStats(String dataset) {
            this.dataset = dataset;
        }
        
        public void calculate() {
            int total = tp + fp + tn + fn;
            coverage = total > 0 ? (double) total / (total + outOfDomain) : 0;
            accuracy = total > 0 ? (double) (tp + tn) / total : 0;
            sensitivity = (tp + fn) > 0 ? (double) tp / (tp + fn) : 0;
            specificity = (tn + fp) > 0 ? (double) tn / (tn + fp) : 0;
            balancedAccuracy = (sensitivity + specificity) / 2;
            ppv = (tp + fp) > 0 ? (double) tp / (tp + fp) : 0;
            npv = (tn + fn) > 0 ? (double) tn / (tn + fn) : 0;
            
            double denominator = Math.sqrt((double)(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
            mcc = denominator > 0 ? (double)(tp * tn - fp * fn) / denominator : 0;
        }
    }




    private final static String TRAINING_SET_PATH = "/data/stero_eds_TS.txt";

    private final PubchemFingerprinter fingerprinter;
    private final SmilesParser smilesParser;
    private final IChemObjectBuilder builder;
    private List<MoleculeData> trainingData;


    public MolecularKNNClassifierFixed() throws Exception {
        builder = DefaultChemObjectBuilder.getInstance();
        fingerprinter = new PubchemFingerprinter(builder);
        smilesParser = new SmilesParser(builder);
        LoadTrainingSet(MolecularKNNClassifierFixed.class.getResource(TRAINING_SET_PATH));
    }


    // ADDED - to load internally the training set from a file as resource
    private void LoadTrainingSet(URL DataSource) throws IOException {

        this.trainingData = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(DataSource.openStream()))) {
            String headerLine = br.readLine();
            if (headerLine == null) {
                throw new IOException("Empty file: " + DataSource.getPath());
            }

            String[] headers = headerLine.split("\t");

            // Find column indices
            int dsIndex = -1, idIndex = -1, smilesIndex = -1, targetIndex = -1;
            for (int i = 0; i < headers.length; i++) {
                String header = headers[i].trim().replace("\"", "");
                if (header.equalsIgnoreCase("DS")) dsIndex = i;
                else if (header.equalsIgnoreCase("ID")) idIndex = i;
                else if (header.equalsIgnoreCase("SMILES")) smilesIndex = i;
                else if (header.equalsIgnoreCase("Target")) targetIndex = i;
            }

            if (dsIndex == -1 || idIndex == -1 || smilesIndex == -1 || targetIndex == -1) {
                throw new IOException("Missing required columns in " + DataSource.getPath() +
                        ". Required: DS, ID, SMILES, Target. Found: " + Arrays.toString(headers));
            }

            String line;
            int lineNumber = 1;

            while ((line = br.readLine()) != null) {
                lineNumber++;
                if (line.trim().isEmpty()) continue;

                String[] parts = line.split("\t");

                if (parts.length > Math.max(Math.max(dsIndex, idIndex), Math.max(smilesIndex, targetIndex))) {
                    try {
                        String ds = parts[dsIndex].trim().replace("\"", "");
                        String id = parts[idIndex].trim().replace("\"", "");
                        String smiles = parts[smilesIndex].trim().replace("\"", "");
                        String target = parts[targetIndex].trim().replace("\"", "");

                        if (smiles.isEmpty()) {
                            System.err.println("Empty SMILES at line " + lineNumber + " in " + DataSource.getPath());
                            continue;
                        }

                        MoleculeData mol = new MoleculeData(ds, id, smiles, target);
                        this.trainingData.add(mol);

                    } catch (ArrayIndexOutOfBoundsException e) {
                        System.err.println("Array index error at line " + lineNumber +
                                " in " + DataSource.getPath() + ": insufficient columns. Expected columns: " +
                                Math.max(Math.max(dsIndex, idIndex), Math.max(smilesIndex, targetIndex)) +
                                ", found: " + parts.length);
                    } catch (Exception e) {
                        System.err.println("Error processing line " + lineNumber +
                                " in " + DataSource.getPath() + ": " + e.getMessage());
                    }
                } else {
                    System.err.println("Insufficient columns at line " + lineNumber + " in " + DataSource.getPath() +
                            ". Expected: " + Math.max(Math.max(dsIndex, idIndex), Math.max(smilesIndex, targetIndex)) +
                            ", found: " + parts.length);
                }
            }
        }

        if (this.trainingData.isEmpty()) {
            throw new IOException("No valid data loaded from " + DataSource.getPath());
        }
    }


    private BitSet calculateFingerprint(String smiles) throws CDKException {
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        
        // Add hydrogens if needed
        try {
            org.openscience.cdk.tools.CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
        } catch (Exception e) {
            // Continue without adding hydrogens if it fails
        }
        
        return fingerprinter.getBitFingerprint(molecule).asBitSet();
    }


    private double calculateSimilarity(BitSet fp1, BitSet fp2) {
        try {
            // Manual calculation of Tanimoto similarity
            BitSet intersection = (BitSet) fp1.clone();
            intersection.and(fp2);
            int intersectionCount = intersection.cardinality();
            
            BitSet union = (BitSet) fp1.clone();
            union.or(fp2);
            int unionCount = union.cardinality();
            
            if (unionCount == 0) return 0.0;
            return (double) intersectionCount / unionCount;
            
        } catch (Exception e) {
            System.err.println("Error calculating similarity: " + e.getMessage());
            return 0.0;
        }
    }


    public Prediction predict(String SMILES) throws Exception {
        MoleculeData testMolecule = new MoleculeData("", "", SMILES, "");
        List<MoleculeData> testData = new ArrayList<>();
        testData.add(testMolecule);
        List<Prediction> res = predict(testData);
        if (res.size()>0)
            return res.get(0);
        throw new Exception("No results returned");
    }


    public List<Prediction> predict(List<MoleculeData> testData) throws Exception {
        List<Prediction> predictions = new ArrayList<>();

        for (MoleculeData testMol : testData) {

            if (testMol.fingerprint == null)
                try {
                    testMol.fingerprint = calculateFingerprint(testMol.smiles);
                } catch (Exception e) {
                    predictions.add(new Prediction(testMol.ds, testMol.id, testMol.smiles,
                            testMol.target, "NA", 0, 0, 0, 0));
                    continue;
                }

            List<Map.Entry<MoleculeData, Double>> similarities = new ArrayList<>();
            
            // Calculate similarities with training data
            for (MoleculeData trainMol : trainingData) {

                // calculate here FP for TS molecule and store it
                if (trainMol.fingerprint == null)
                    try {
                        trainMol.fingerprint = calculateFingerprint(trainMol.smiles);
                    } catch (Exception e) {
                        throw new Exception("Unable to calculate FP for training set molecule with SMILES " + trainMol.smiles);
                    }

                double similarity = calculateSimilarity(testMol.fingerprint, trainMol.fingerprint);
                similarities.add(new AbstractMap.SimpleEntry<>(trainMol, similarity));
            }
            
            // Sort by similarity (descending)
            similarities.sort((a, b) -> Double.compare(b.getValue(), a.getValue()));
            
            // Filter by minimum similarity and take top k
            List<Map.Entry<MoleculeData, Double>> topK = similarities.stream()
                .filter(entry -> entry.getValue() >= MIN_SIMILARITY)
                .limit(K)
                .collect(Collectors.toList());
            
            String prediction = "OutOfDomain";
            double avgSimilarity = 0.0;
            int neighbors = topK.size();
            int positives = 0;
            int negatives = 0;
            
            if (!topK.isEmpty()) {
                // Count votes and positive/negative neighbors
                Map<String, Integer> votes = new HashMap<>();
                double totalSimilarity = 0.0;
                
                for (Map.Entry<MoleculeData, Double> entry : topK) {
                    String target = entry.getKey().target;
                    votes.put(target, votes.getOrDefault(target, 0) + 1);
                    totalSimilarity += entry.getValue();
                    
                    // Count positives and negatives
                    boolean isPositive = target.equals("1") || target.equalsIgnoreCase("positive") || target.equalsIgnoreCase("active");
                    if (isPositive) {
                        positives++;
                    } else {
                        negatives++;
                    }
                }
                
                // Get majority vote
                prediction = votes.entrySet().stream()
                    .max(Map.Entry.comparingByValue())
                    .map(Map.Entry::getKey)
                    .orElse("OutOfDomain");
                
                avgSimilarity = totalSimilarity / topK.size();
            }
            
            predictions.add(new Prediction(testMol.ds, testMol.id, testMol.smiles, 
                                         testMol.target, prediction, avgSimilarity, neighbors, positives, negatives));
        }
        
        return predictions;
    }


    private ClassificationStats calculateStats(List<Prediction> predictions, String datasetName) {
        ClassificationStats stats = new ClassificationStats(datasetName);
        
        for (Prediction pred : predictions) {
            if (pred.prediction.equals("OutOfDomain")) {
                stats.outOfDomain++;
            } else {
                boolean actualPositive = pred.target.equals("1") || pred.target.equalsIgnoreCase("positive") || pred.target.equalsIgnoreCase("active");
                boolean predictedPositive = pred.prediction.equals("1") || pred.prediction.equalsIgnoreCase("positive") || pred.prediction.equalsIgnoreCase("active");
                
                if (actualPositive && predictedPositive) stats.tp++;
                else if (!actualPositive && predictedPositive) stats.fp++;
                else if (!actualPositive && !predictedPositive) stats.tn++;
                else stats.fn++;
            }
        }
        
        stats.calculate();
        return stats;
    }



}