package insilico.dio1.model.modelKNN;

import insilico.core.exception.InitFailureException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.stream.Collectors;

public class DioMolecularKNNClassifierFixed {
    
    // Fixed parameters
    private static final int K = 5;
    private static final double MIN_SIMILARITY = 0.8;
    
    private static class MoleculeData {
        String ds, id, smiles, target;
        BitSet fingerprint;
        
        public MoleculeData(String ds, String id, String smiles, String target) {
            this.ds = ds;
            this.id = id;
            this.smiles = smiles;
            this.target = target;
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
    
    private PubchemFingerprinter fingerprinter;
    private SmilesParser smilesParser;
    private List<MoleculeData> trainingData;
    private IChemObjectBuilder builder;
    private final static String TS_URL = "/data/dio1_TS.txt";
    
    public DioMolecularKNNClassifierFixed() throws InitFailureException {
        try {
            builder = DefaultChemObjectBuilder.getInstance();
            fingerprinter = new PubchemFingerprinter(builder);
            smilesParser = new SmilesParser(builder);

            URL DataSource = DioMolecularKNNClassifierFixed.class.getResource(TS_URL);
            this.trainingData = loadData(DataSource);
        } catch (Throwable e) {
            throw new InitFailureException(e);
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
    
    private List<MoleculeData> loadData(URL source) throws IOException {
        List<MoleculeData> data = new ArrayList<>();

        String filename = source.getPath();


        try (BufferedReader br = new BufferedReader(new InputStreamReader(source.openStream()))) {
            String headerLine = br.readLine();
            if (headerLine == null) {
                throw new IOException("Empty file as input" + filename);
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
                throw new IOException("Missing required columns in " + filename + 
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
                        String smiles = parts[smilesIndex].trim().replace("\"", "");
                        if (smiles.isEmpty()) {
                            System.err.println("Empty SMILES at line " + lineNumber + " in " + filename);
                            continue;
                        }
                        
                        MoleculeData mol = new MoleculeData(
                            parts[dsIndex].trim().replace("\"", ""),
                            parts[idIndex].trim().replace("\"", ""),
                            smiles,
                            parts[targetIndex].trim().replace("\"", "")
                        );
                        
                        try {
                            mol.fingerprint = calculateFingerprint(mol.smiles);
                            data.add(mol);
                        } catch (CDKException e) {
                            System.err.println("Error calculating fingerprint for SMILES at line " + lineNumber + 
                                " in " + filename + ": " + smiles + " - " + e.getMessage());
                            continue;
                        }
                    } catch (ArrayIndexOutOfBoundsException e) {
                        System.err.println("Array index error at line " + lineNumber + 
                            " in " + filename + ": insufficient columns");
                    } catch (Exception e) {
                        System.err.println("Error processing line " + lineNumber + 
                            " in " + filename + ": " + e.getMessage());
                    }
                } else {
                    System.err.println("Insufficient columns at line " + lineNumber + " in " + filename);
                }
            }
        }
        
        if (data.isEmpty()) {
            throw new IOException("No valid data loaded from " + filename);
        }
        
        return data;
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

    public Prediction predict(String SMILES) {

        MoleculeData mol = new MoleculeData("", "", SMILES, "");
        try {
            mol.fingerprint = calculateFingerprint(SMILES);

            ArrayList<MoleculeData> bufInput = new ArrayList<>();
            bufInput.add(mol);
            List<Prediction> predictions = predict(bufInput);
            return predictions.get(0);

        } catch (CDKException e) {
            return new Prediction("", "", "", "", "-999", 0, 0, 0, 0);
        }

    }


    private List<Prediction> predict(List<MoleculeData> testData) {
        List<Prediction> predictions = new ArrayList<>();
        int processedCount = 0;
        
        for (MoleculeData testMol : testData) {
            processedCount++;
            if (processedCount % 100 == 0) {
                System.out.printf("Processed %d/%d molecules%n", processedCount, testData.size());
            }
            
            List<Map.Entry<MoleculeData, Double>> similarities = new ArrayList<>();
            
            // Calculate similarities with training data
            for (MoleculeData trainMol : trainingData) {
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
    
//    public void run() throws IOException {
//        // Check if files exist
//        if (!new File("TS.csv").exists()) {
//            throw new FileNotFoundException("File TS.csv not found in current directory");
//        }
//        if (!new File("CS.csv").exists()) {
//            throw new FileNotFoundException("File CS.csv not found in current directory");
//        }
//        if (!new File("VS.csv").exists()) {
//            throw new FileNotFoundException("File VS.csv not found in current directory");
//        }
//
//        System.out.printf("Using fixed parameters: K=%d, Min Similarity=%.1f%n", K, MIN_SIMILARITY);
//
//        // Load data
//        System.out.println("Loading training data...");
//        trainingData = loadData("TS.csv");
//        System.out.println("Loaded " + trainingData.size() + " training molecules");
//
//        System.out.println("Loading calibration data...");
//        List<MoleculeData> calibrationData = loadData("CS.csv");
//        System.out.println("Loaded " + calibrationData.size() + " calibration molecules");
//
//        System.out.println("Loading validation data...");
//        List<MoleculeData> validationData = loadData("VS.csv");
//        System.out.println("Loaded " + validationData.size() + " validation molecules");
//
//        // Generate predictions with fixed parameters
//        System.out.println("Generating predictions for training set...");
//        List<Prediction> tsPredictions = predict(trainingData);
//
//        System.out.println("Generating predictions for calibration set...");
//        List<Prediction> csPredictions = predict(calibrationData);
//
//        System.out.println("Generating predictions for validation set...");
//        List<Prediction> vsPredictions = predict(validationData);
//
//        // Combine all predictions
//        List<Prediction> allPredictions = new ArrayList<>();
//        allPredictions.addAll(tsPredictions);
//        allPredictions.addAll(csPredictions);
//        allPredictions.addAll(vsPredictions);
//
//        System.out.println("Saving predictions...");
//        savePredictions(allPredictions);
//
//        // Calculate and save statistics
//        System.out.println("Calculating statistics...");
//        ClassificationStats tsStats = calculateStats(tsPredictions, "TS");
//        ClassificationStats csStats = calculateStats(csPredictions, "CS");
//        ClassificationStats vsStats = calculateStats(vsPredictions, "VS");
//
//        List<ClassificationStats> finalStats = Arrays.asList(tsStats, csStats, vsStats);
//        saveStats(finalStats);
//
//        // Print summary
//        System.out.println("\n=== RESULTS SUMMARY ===");
//        System.out.printf("Parameters used: K=%d, Min Similarity=%.1f%n", K, MIN_SIMILARITY);
//        System.out.println("\nDataset Statistics:");
//        for (ClassificationStats stat : finalStats) {
//            System.out.printf("%s: Coverage=%.3f, MCC=%.3f, Accuracy=%.3f, Sensitivity=%.3f, Specificity=%.3f%n",
//                stat.dataset, stat.coverage, stat.mcc, stat.accuracy, stat.sensitivity, stat.specificity);
//        }
//
//        System.out.println("\nFiles generated:");
//        System.out.println("- predictions_fixed.csv: All predictions");
//        System.out.println("- stats_fixed.csv: Detailed statistics");
//        System.out.println("\nAnalysis complete!");
//    }
    
    private void savePredictions(List<Prediction> predictions) throws IOException {
        // Set locale to ensure decimal point separator
        Locale.setDefault(Locale.US);
        
        try (PrintWriter pw = new PrintWriter(new FileWriter("predictions_fixed.csv"))) {
            pw.println("DS,ID,SMILES,Target,Prediction,Similarity,Neighbors,Positives,Negatives");
            for (Prediction pred : predictions) {
                pw.printf(Locale.US, "%s,%s,%s,%s,%s,%.4f,%d,%d,%d%n",
                         pred.ds, pred.id, pred.smiles, pred.target, pred.prediction, 
                         pred.similarity, pred.neighbors, pred.positives, pred.negatives);
            }
        }
    }
    
    private void saveStats(List<ClassificationStats> stats) throws IOException {
        // Set locale to ensure decimal point separator
        Locale.setDefault(Locale.US);
        
        try (PrintWriter pw = new PrintWriter(new FileWriter("stats_fixed.csv"))) {
            pw.println("Dataset,TP,FP,TN,FN,OutOfDomain,Coverage,BalancedAccuracy,Accuracy,Sensitivity,Specificity,PPV,NPV,MCC");
            for (ClassificationStats stat : stats) {
                pw.printf(Locale.US, "%s,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f%n",
                         stat.dataset, stat.tp, stat.fp, stat.tn, stat.fn, stat.outOfDomain,
                         stat.coverage, stat.balancedAccuracy, stat.accuracy, stat.sensitivity,
                         stat.specificity, stat.ppv, stat.npv, stat.mcc);
            }
        }
    }
    
//    public static void main(String[] args) {
//        try {
//            System.out.println("Starting Molecular kNN Classifier (Fixed Parameters)...");
//            System.out.println("Current working directory: " + System.getProperty("user.dir"));
//
//            MolecularKNNClassifierFixed classifier = new MolecularKNNClassifierFixed();
//            classifier.run();
//
//            System.out.println("Analysis completed successfully!");
//        } catch (FileNotFoundException e) {
//            System.err.println("File not found: " + e.getMessage());
//            System.err.println("Please ensure TS.csv, CS.csv, and VS.csv are in the current directory");
//        } catch (IOException e) {
//            System.err.println("IO Error: " + e.getMessage());
//            e.printStackTrace();
//        } catch (Exception e) {
//            System.err.println("Unexpected error: " + e.getMessage());
//            e.printStackTrace();
//        }
//    }
}