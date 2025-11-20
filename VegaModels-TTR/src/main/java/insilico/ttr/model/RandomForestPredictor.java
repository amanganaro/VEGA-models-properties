package insilico.ttr.model;

import insilico.core.exception.InitFailureException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RandomForestPredictor {
    
    // Classe per rappresentare un nodo dell'albero decisionale
    static class TreeNode {
        String id;
        String score;
        String field;
        String operator;
        double value;
        List<TreeNode> children;
        Map<String, Double> scoreDistribution;
        boolean isLeaf;
        
        TreeNode() {
            children = new ArrayList<>();
            scoreDistribution = new HashMap<>();
            isLeaf = false;
        }
    }
    
    // Classe per rappresentare un record di dati
    static class DataRecord {
        String ds;
        String id;
        String smiles;
        String target;
        double mats1e;
        double ssoh;
        
        DataRecord(String ds, String id, String smiles, String target, double mats1e, double ssoh) {
            this.ds = ds;
            this.id = id;
            this.smiles = smiles;
            this.target = target;
            this.mats1e = mats1e;
            this.ssoh = ssoh;
        }
    }
    
    // Classe per rappresentare una predizione
    public static class Prediction {
        DataRecord record;
        public String prediction;
        public double confidence;
        
        Prediction(DataRecord record, String prediction, double confidence) {
            this.record = record;
            this.prediction = prediction;
            this.confidence = confidence;
        }
    }
    
    // Classe per le statistiche
    static class Statistics {
        int tp, fp, tn, fn, outOfDomain;
        double balancedAccuracy, accuracy, sensitivity, specificity, ppv, npv, mcc;
        
        void calculate() {
            int total = tp + fp + tn + fn;
            if (total == 0) return;
            
            accuracy = (double)(tp + tn) / total;
            sensitivity = tp + fn > 0 ? (double)tp / (tp + fn) : 0.0;
            specificity = tn + fp > 0 ? (double)tn / (tn + fp) : 0.0;
            ppv = tp + fp > 0 ? (double)tp / (tp + fp) : 0.0;
            npv = tn + fn > 0 ? (double)tn / (tn + fn) : 0.0;
            balancedAccuracy = (sensitivity + specificity) / 2.0;
            
            double denominator = Math.sqrt((double)(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
            mcc = denominator > 0 ? (double)(tp * tn - fp * fn) / denominator : 0.0;
        }
    }
    
    private final List<TreeNode> trees;
    private final TtrDescriptors descriptors;


    private final static String PMML_PATH = "/data/ttr_rf_model.pmml";

    public RandomForestPredictor() throws InitFailureException {
        trees = new ArrayList<>();
        descriptors = new TtrDescriptors();
        try {
            loadModel(RandomForestPredictor.class.getResource(PMML_PATH));
        } catch (Exception e) {
            throw new InitFailureException("unable to load PMML file - " + e.getMessage());
        }
    }


    // Carica il modello PMML
    private void loadModel(URL pmmlFile) throws Exception {

        InputStream is = pmmlFile.openStream();

        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
//        Document doc = builder.parse(new File(pmmlFile));
        Document doc = builder.parse(is);

        NodeList segments = doc.getElementsByTagName("Segment");
        
        for (int i = 0; i < segments.getLength(); i++) {
            Element segment = (Element) segments.item(i);
            NodeList treeModels = segment.getElementsByTagName("TreeModel");
            
            if (treeModels.getLength() > 0) {
                Element treeModel = (Element) treeModels.item(0);
                NodeList nodes = treeModel.getElementsByTagName("Node");
                
                if (nodes.getLength() > 0) {
                    TreeNode rootNode = parseNode((Element) nodes.item(0), nodes);
                    trees.add(rootNode);
                }
            }
        }
        
    }
    
    // Parsing ricorsivo dei nodi dell'albero
    private TreeNode parseNode(Element nodeElement, NodeList allNodes) {
        TreeNode node = new TreeNode();
        node.id = nodeElement.getAttribute("id");
        node.score = nodeElement.getAttribute("score");
        
        // Parse delle distribuzioni di score
        NodeList scoreDistributions = nodeElement.getElementsByTagName("ScoreDistribution");
        for (int i = 0; i < scoreDistributions.getLength(); i++) {
            Element scoreDist = (Element) scoreDistributions.item(i);
            String value = scoreDist.getAttribute("value");
            double recordCount = Double.parseDouble(scoreDist.getAttribute("recordCount"));
            node.scoreDistribution.put(value, recordCount);
        }
        
        // Parse dei predicati
        NodeList simplePredicates = nodeElement.getElementsByTagName("SimplePredicate");
        if (simplePredicates.getLength() > 0) {
            Element predicate = (Element) simplePredicates.item(0);
            node.field = predicate.getAttribute("field");
            node.operator = predicate.getAttribute("operator");
            String valueStr = predicate.getAttribute("value");
            if (!valueStr.isEmpty()) {
                node.value = Double.parseDouble(valueStr);
            }
        }
        
        // Trova i nodi figli
        List<String> childIds = new ArrayList<>();
        NodeList childNodes = nodeElement.getChildNodes();
        for (int i = 0; i < childNodes.getLength(); i++) {
            if (childNodes.item(i).getNodeType() == Node.ELEMENT_NODE) {
                Element child = (Element) childNodes.item(i);
                if ("Node".equals(child.getTagName())) {
                    childIds.add(child.getAttribute("id"));
                }
            }
        }
        
        // Parse ricorsivo dei figli
        for (String childId : childIds) {
            for (int i = 0; i < allNodes.getLength(); i++) {
                Element childElement = (Element) allNodes.item(i);
                if (childId.equals(childElement.getAttribute("id"))) {
                    TreeNode childNode = parseNode(childElement, allNodes);
                    node.children.add(childNode);
                    break;
                }
            }
        }
        
        node.isLeaf = node.children.isEmpty();
        return node;
    }
    
    // Predizione per un singolo record usando un albero
    private String predictWithTree(TreeNode root, DataRecord record) {
        TreeNode current = root;
        
        while (!current.isLeaf && !current.children.isEmpty()) {
            boolean conditionMet = false;
            
            for (TreeNode child : current.children) {
                if (evaluatePredicate(child, record)) {
                    current = child;
                    conditionMet = true;
                    break;
                }
            }
            
            if (!conditionMet) {
                break;
            }
        }
        
        return current.score;
    }
    
    // Valutazione del predicato
    private boolean evaluatePredicate(TreeNode node, DataRecord record) {
        if (node.field == null) return true;
        
        double fieldValue;
        switch (node.field) {
            case "MATS1e":
                fieldValue = record.mats1e;
                break;
            case "SsOH":
                fieldValue = record.ssoh;
                break;
            default:
                return true;
        }
        
        switch (node.operator) {
            case "lessOrEqual":
                return fieldValue <= node.value;
            case "greaterThan":
                return fieldValue > node.value;
            case "isMissing":
                return Double.isNaN(fieldValue);
            default:
                return true;
        }
    }

    public Prediction predict(String SMILES) throws Exception {
        double[] desc = descriptors.Scale(descriptors.Calculate(SMILES));
        DataRecord dr = new DataRecord("", "", SMILES, "", desc[0], desc[1]);
        Prediction pred = predict(dr);
        return pred;
    }

    // Predizione con Random Forest
    public Prediction predict(DataRecord record) {
        Map<String, Integer> votes = new HashMap<>();
        
        for (TreeNode tree : trees) {
            String prediction = predictWithTree(tree, record);
            votes.put(prediction, votes.getOrDefault(prediction, 0) + 1);
        }
        
        String finalPrediction = "";
        int maxVotes = 0;
        for (Map.Entry<String, Integer> entry : votes.entrySet()) {
            if (entry.getValue() > maxVotes) {
                maxVotes = entry.getValue();
                finalPrediction = entry.getKey();
            }
        }
        
        double confidence = (double) maxVotes / trees.size();
        return new Prediction(record, finalPrediction, confidence);
    }



//
//    // Caricamento dei dati da CSV
//    public List<DataRecord> loadCSV(String filename) throws IOException {
//        List<DataRecord> records = new ArrayList<>();
//
//        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
//            String line = br.readLine(); // Skip header
//
//            while ((line = br.readLine()) != null) {
//                String[] parts = line.split(",");
//                if (parts.length >= 6) {
//                    String ds = parts[0].replace("\"", "");
//                    String id = parts[1].replace("\"", "");
//                    String smiles = parts[2].replace("\"", "");
//                    String target = parts[3].replace("\"", "");
//                    double mats1e = Double.parseDouble(parts[4]);
//                    double ssoh = Double.parseDouble(parts[5]);
//
//                    records.add(new DataRecord(ds, id, smiles, target, mats1e, ssoh));
//                }
//            }
//        }
//
//        return records;
//    }
//
//    // Salvataggio delle predizioni
//    public void savePredictions(List<Prediction> predictions, String filename) throws IOException {
//        try (PrintWriter pw = new PrintWriter(new FileWriter(filename))) {
//            pw.println("DS,ID,SMILES,Target,MATS1e,SsOH,Prediction,Confidence");
//
//            for (Prediction pred : predictions) {
//                DataRecord r = pred.record;
//                pw.printf("%s,%s,\"%s\",%s,%.6f,%.6f,%s,%.6f%n",
//                    r.ds, r.id, r.smiles, r.target, r.mats1e, r.ssoh,
//                    pred.prediction, pred.confidence);
//            }
//        }
//    }
//
//    // Calcolo delle statistiche
//    public Map<String, Statistics> calculateStatistics(List<Prediction> predictions) {
//        Map<String, Statistics> stats = new HashMap<>();
//
//        // Raggruppa per dataset
//        Map<String, List<Prediction>> byDataset = new HashMap<>();
//        for (Prediction pred : predictions) {
//            String ds = pred.record.ds;
//            byDataset.computeIfAbsent(ds, k -> new ArrayList<>()).add(pred);
//        }
//
//        // Calcola statistiche per ogni dataset
//        for (Map.Entry<String, List<Prediction>> entry : byDataset.entrySet()) {
//            Statistics stat = new Statistics();
//
//            for (Prediction pred : entry.getValue()) {
//                String actual = pred.record.target;
//                String predicted = pred.prediction;
//
//                if ("1".equals(actual) && "1".equals(predicted)) {
//                    stat.tp++;
//                } else if ("0".equals(actual) && "1".equals(predicted)) {
//                    stat.fp++;
//                } else if ("0".equals(actual) && "0".equals(predicted)) {
//                    stat.tn++;
//                } else if ("1".equals(actual) && "0".equals(predicted)) {
//                    stat.fn++;
//                } else {
//                    stat.outOfDomain++;
//                }
//            }
//
//            stat.calculate();
//            stats.put(entry.getKey(), stat);
//        }
//
//        return stats;
//    }
//
//    // Salvataggio delle statistiche
//    public void saveStatistics(Map<String, Statistics> stats, String filename) throws IOException {
//        try (PrintWriter pw = new PrintWriter(new FileWriter(filename))) {
//            pw.println("Dataset,TP,FP,TN,FN,OutOfDomain,Accuracy,BalancedAccuracy,Sensitivity,Specificity,PPV,NPV,MCC");
//
//            for (Map.Entry<String, Statistics> entry : stats.entrySet()) {
//                Statistics s = entry.getValue();
//                pw.printf("%s,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f%n",
//                    entry.getKey(), s.tp, s.fp, s.tn, s.fn, s.outOfDomain,
//                    s.accuracy, s.balancedAccuracy, s.sensitivity, s.specificity,
//                    s.ppv, s.npv, s.mcc);
//            }
//        }
//    }
    
//    public static void main(String[] args) {
//        try {
//            RandomForestPredictor predictor = new RandomForestPredictor();
//
//            // Carica il modello PMML
//            predictor.loadModel("model.pmml");
//
//            // Lista per tutte le predizioni
//            List<Prediction> allPredictions = new ArrayList<>();
//
//            // Carica e predici sui diversi dataset
//            String[] datasetFiles = {"TS.csv", "CS.csv", "VS.csv"};
//
//            for (String file : datasetFiles) {
//                File f = new File(file);
//                if (f.exists()) {
//                    System.out.println("Elaborando " + file + "...");
//                    List<DataRecord> records = predictor.loadCSV(file);
//
//                    for (DataRecord record : records) {
//                        Prediction pred = predictor.predict(record);
//                        allPredictions.add(pred);
//                    }
//
//                    System.out.println("Elaborate " + records.size() + " righe da " + file);
//                } else {
//                    System.out.println("File " + file + " non trovato, saltato.");
//                }
//            }
//
//            // Salva le predizioni
//            predictor.savePredictions(allPredictions, "predictions.csv");
//            System.out.println("Predizioni salvate in predictions.csv");
//
//            // Calcola e salva le statistiche
//            Map<String, Statistics> statistics = predictor.calculateStatistics(allPredictions);
//            predictor.saveStatistics(statistics, "statistics.csv");
//            System.out.println("Statistiche salvate in statistics.csv");
//
//            // Stampa riepilogo
//            System.out.println("\nRiepilogo:");
//            for (Map.Entry<String, Statistics> entry : statistics.entrySet()) {
//                Statistics s = entry.getValue();
//                System.out.printf("Dataset %s: Accuracy=%.3f, Balanced Accuracy=%.3f, MCC=%.3f%n",
//                    entry.getKey(), s.accuracy, s.balancedAccuracy, s.mcc);
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
}
