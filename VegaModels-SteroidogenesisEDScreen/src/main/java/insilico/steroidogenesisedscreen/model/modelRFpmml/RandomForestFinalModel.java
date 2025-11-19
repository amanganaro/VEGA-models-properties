package insilico.steroidogenesisedscreen.model.modelRFpmml;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/*
  RandomForestFinalModel.java con esportazione PMML
  - Esporta il modello in formato PMML (XML standard per modelli predittivi)
  - Mantiene la serializzazione Java (.ser) come backup
  - PMML è facilmente importabile in altri linguaggi (Python, R, etc.)
*/

/* -------------------------
   Classifier (Random Forest)
   ------------------------- */
class RandomForestClassifier implements Serializable {
    private static final long serialVersionUID = 1L;

    // Package-private per permettere l'accesso da PMMLExporter
    List<DecisionTree> trees;
    String[] featureNames;
    Map<String, Integer> classToIndex;
    String[] indexToClass;
    int numTrees;
    int maxDepth;
    int minNodeSize;
    Random random;

    public RandomForestClassifier(int numTrees, int maxDepth, int minNodeSize) {
        this.numTrees = numTrees;
        this.maxDepth = maxDepth;
        this.minNodeSize = minNodeSize;
        this.trees = new ArrayList<>();
        this.random = new Random(17);
    }

    public void train(DataSet trainingData) {
        this.featureNames = trainingData.getFeatureNames();
        this.classToIndex = trainingData.getClassMapping();
        this.indexToClass = new String[classToIndex.size()];
        for (Map.Entry<String, Integer> entry : classToIndex.entrySet()) {
            indexToClass[entry.getValue()] = entry.getKey();
        }

        List<Instance> instances = new ArrayList<>(trainingData.getInstances());

        System.out.println("Costruzione Random Forest con " + numTrees + " alberi...");
        for (int t = 0; t < numTrees; t++) {
            List<Instance> sample = new ArrayList<>();
            for (int i = 0; i < instances.size(); i++) {
                int idx = random.nextInt(instances.size());
                sample.add(instances.get(idx));
            }

            DecisionTree tree = new DecisionTree(maxDepth, minNodeSize, random);
            tree.buildTree(sample, featureNames);
            trees.add(tree);
            if ((t + 1) % 50 == 0) {
                System.out.println("  Alberi costruiti: " + (t + 1) + "/" + numTrees);
            }
        }

        System.out.println("Training completato.");
    }

    public PredictionResult predict(Instance instance) {
        Map<String, Integer> votes = new HashMap<>();
        Map<String, Double> confidences = new HashMap<>();

        for (DecisionTree tree : trees) {
            String pred = tree.predict(instance);
            votes.put(pred, votes.getOrDefault(pred, 0) + 1);
        }

        String bestClass = null;
        int bestVotes = -1;
        for (Map.Entry<String, Integer> e : votes.entrySet()) {
            if (e.getValue() > bestVotes) {
                bestVotes = e.getValue();
                bestClass = e.getKey();
            }
        }

        double confidence = (double) bestVotes / (double) trees.size();
        return new PredictionResult(bestClass, confidence);
    }

    public List<PredictionResult> predict(DataSet dataSet) {
        List<PredictionResult> predictions = new ArrayList<>();
        System.out.println("Predizione su " + dataSet.getInstances().size() + " istanze...");

        for (int i = 0; i < dataSet.getInstances().size(); i++) {
            predictions.add(predict(dataSet.getInstances().get(i)));

            if ((i + 1) % 1000 == 0) {
                System.out.println("Predette " + (i + 1) + "/" + dataSet.getInstances().size() + " istanze");
            }
        }

        return predictions;
    }

    // Getters per esportazione PMML
    public int getNumTrees() { return numTrees; }
    public List<DecisionTree> getTrees() { return trees; }
    public String[] getFeatureNames() { return featureNames; }
    public String[] getIndexToClass() { return indexToClass; }
}

/* -------------------------
   Decision Tree
   ------------------------- */
class DecisionTree implements Serializable {
    private static final long serialVersionUID = 1L;

    // Package-private per permettere l'accesso da PMMLExporter
    TreeNode root;
    int maxDepth;
    int minNodeSize;
    Random random;

    public DecisionTree(int maxDepth, int minNodeSize, Random random) {
        this.maxDepth = maxDepth;
        this.minNodeSize = minNodeSize;
        this.random = random;
    }

    public void buildTree(List<Instance> instances, String[] featureNames) {
        this.root = buildTreeRecursive(instances, 0, featureNames);
    }

    private TreeNode buildTreeRecursive(List<Instance> instances, int depth, String[] featureNames) {
        if (instances.isEmpty()) return null;

        Map<String, Integer> classCounts = new HashMap<>();
        for (Instance instance : instances) {
            classCounts.put(instance.getTarget(),
                classCounts.getOrDefault(instance.getTarget(), 0) + 1);
        }

        String majorityClass = classCounts.entrySet().stream()
            .max(Map.Entry.comparingByValue())
            .get().getKey();

        if (classCounts.size() == 1 || depth >= maxDepth || instances.size() < minNodeSize) {
            return new TreeNode(majorityClass);
        }

        Split bestSplit = findBestSplit(instances, featureNames);
        if (bestSplit == null) {
            return new TreeNode(majorityClass);
        }

        List<Instance> left = new ArrayList<>();
        List<Instance> right = new ArrayList<>();
        for (Instance inst : instances) {
            if (inst.getFeatureValue(bestSplit.featureIndex) <= bestSplit.threshold) {
                left.add(inst);
            } else {
                right.add(inst);
            }
        }

        TreeNode leftChild = buildTreeRecursive(left, depth + 1, featureNames);
        TreeNode rightChild = buildTreeRecursive(right, depth + 1, featureNames);
        return new TreeNode(bestSplit.featureIndex, bestSplit.threshold, leftChild, rightChild);
    }

    private Split findBestSplit(List<Instance> instances, String[] featureNames) {
        Split bestSplit = null;
        double bestGini = Double.MAX_VALUE;

        int numFeaturesToTry = Math.max(1, (int) Math.sqrt(featureNames.length));
        List<Integer> featureIndices = new ArrayList<>();
        for (int i = 0; i < featureNames.length; i++) {
            featureIndices.add(i);
        }
        Collections.shuffle(featureIndices, random);

        for (int fi = 0; fi < numFeaturesToTry; fi++) {
            int featureIndex = featureIndices.get(fi);

            Set<Double> values = new HashSet<>();
            for (Instance inst : instances) {
                values.add(inst.getFeatureValue(featureIndex));
            }

            List<Double> sorted = new ArrayList<>(values);
            Collections.sort(sorted);

            for (int i = 0; i < sorted.size() - 1; i++) {
                double threshold = (sorted.get(i) + sorted.get(i + 1)) / 2.0;
                double gini = computeWeightedGini(instances, featureIndex, threshold);
                if (gini < bestGini) {
                    bestGini = gini;
                    bestSplit = new Split(featureIndex, threshold);
                }
            }
        }

        return bestSplit;
    }

    private double computeWeightedGini(List<Instance> instances, int featureIndex, double threshold) {
        List<Instance> left = new ArrayList<>();
        List<Instance> right = new ArrayList<>();
        for (Instance inst : instances) {
            if (inst.getFeatureValue(featureIndex) <= threshold) {
                left.add(inst);
            } else {
                right.add(inst);
            }
        }

        if (left.isEmpty() || right.isEmpty()) return Double.MAX_VALUE;

        double totalSize = instances.size();
        double leftWeight = left.size() / totalSize;
        double rightWeight = right.size() / totalSize;

        return leftWeight * calculateGini(left) + rightWeight * calculateGini(right);
    }

    private double calculateGini(List<Instance> instances) {
        if (instances.isEmpty()) return 0.0;

        Map<String, Integer> classCounts = new HashMap<>();
        for (Instance instance : instances) {
            classCounts.put(instance.getTarget(),
                classCounts.getOrDefault(instance.getTarget(), 0) + 1);
        }

        double gini = 1.0;
        double total = instances.size();

        for (int count : classCounts.values()) {
            double probability = count / total;
            gini -= probability * probability;
        }

        return gini;
    }

    public String predict(Instance instance) {
        TreeNode node = root;
        while (node != null && !node.isLeaf()) {
            double v = instance.getFeatureValue(node.getFeatureIndex());
            if (v <= node.getThreshold()) {
                node = node.getLeft();
            } else {
                node = node.getRight();
            }
        }
        return node == null ? null : node.getPrediction();
    }

    public TreeNode getRoot() { return root; }
    
    // Setter per permettere l'import da PMML
    public void setRoot(TreeNode root) { this.root = root; }
}

/* -------------------------
   TreeNode, Split
   ------------------------- */
class TreeNode implements Serializable {
    private static final long serialVersionUID = 1L;

    private int featureIndex;
    private double threshold;
    private TreeNode left;
    private TreeNode right;
    private String prediction;

    public TreeNode(String prediction) {
        this.prediction = prediction;
    }

    public TreeNode(int featureIndex, double threshold, TreeNode left, TreeNode right) {
        this.featureIndex = featureIndex;
        this.threshold = threshold;
        this.left = left;
        this.right = right;
    }

    public boolean isLeaf() { return prediction != null; }
    public int getFeatureIndex() { return featureIndex; }
    public double getThreshold() { return threshold; }
    public TreeNode getLeft() { return left; }
    public TreeNode getRight() { return right; }
    public String getPrediction() { return prediction; }
}

class Split implements Serializable {
    private static final long serialVersionUID = 1L;
    public int featureIndex;
    public double threshold;
    public Split(int featureIndex, double threshold) {
        this.featureIndex = featureIndex;
        this.threshold = threshold;
    }
}

/* -------------------------
   TreeStructure (helper per import PMML)
   ------------------------- */
class TreeStructure {
    public TreeNode root;
    
    public TreeStructure() {
        this.root = null;
    }
}

/* -------------------------
   PMML Exporter & Importer
   ------------------------- */
class PMMLExporter {
    
    /**
     * Carica un modello Random Forest da file PMML
     */
    public static RandomForestClassifier importFromPMML(String filename) throws IOException {
        System.out.println("Caricamento modello da PMML: " + filename);
        
        // Leggi il file PMML
        StringBuilder xmlContent = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                xmlContent.append(line).append("\n");
            }
        }
        
        String xml = xmlContent.toString();
        
        // Parsing manuale del PMML (semplice parser XML)
        List<String> featureNames = new ArrayList<>();
        List<String> classNames = new ArrayList<>();
        List<TreeStructure> treeStructures = new ArrayList<>();
        
        // Estrai feature names
        int dataFieldStart = 0;
        while ((dataFieldStart = xml.indexOf("<DataField name=\"", dataFieldStart)) != -1) {
            dataFieldStart += 17; // lunghezza di "<DataField name=\""
            int dataFieldEnd = xml.indexOf("\"", dataFieldStart);
            String fieldName = xml.substring(dataFieldStart, dataFieldEnd);
            
            // Controlla se è il target
            if (xml.indexOf("<Value value=", dataFieldEnd) != -1 && 
                xml.indexOf("<Value value=", dataFieldEnd) < xml.indexOf("</DataField>", dataFieldEnd)) {
                // È il target, estrai i valori delle classi
                int valueStart = dataFieldEnd;
                while ((valueStart = xml.indexOf("<Value value=\"", valueStart)) != -1 && 
                       valueStart < xml.indexOf("</DataField>", dataFieldEnd)) {
                    valueStart += 14;
                    int valueEnd = xml.indexOf("\"", valueStart);
                    if (valueEnd > dataFieldEnd && valueEnd < xml.indexOf("</DataField>", dataFieldEnd)) {
                        classNames.add(xml.substring(valueStart, valueEnd));
                        valueStart = valueEnd;
                    } else {
                        break;
                    }
                }
            } else {
                boolean cond = !fieldName.equals("Target");
                if (cond)
                    featureNames.add(fieldName);
            }
            
            dataFieldStart = dataFieldEnd;
        }
        
        System.out.println("Features estratte: " + featureNames.size());
        System.out.println("Classi estratte: " + classNames);
        
        // Estrai gli alberi
        int segmentStart = 0;
        while ((segmentStart = xml.indexOf("<Segment id=", segmentStart)) != -1) {
            int treeModelStart = xml.indexOf("<TreeModel", segmentStart);
            int treeModelEnd = xml.indexOf("</TreeModel>", treeModelStart);
            
            if (treeModelStart == -1 || treeModelEnd == -1) break;
            
            String treeXml = xml.substring(treeModelStart, treeModelEnd + 12);
            TreeStructure tree = parseTree(treeXml, featureNames);
            treeStructures.add(tree);
            
            segmentStart = treeModelEnd;
        }
        
        System.out.println("Alberi estratti: " + treeStructures.size());
        
        // Ricostruisci il modello
        RandomForestClassifier model = new RandomForestClassifier(
            treeStructures.size(), 
            100, // maxDepth (non critico per modello già addestrato)
            1    // minNodeSize (non critico per modello già addestrato)
        );
        
        // Imposta i parametri del modello
        model.featureNames = featureNames.toArray(new String[0]);
        model.indexToClass = classNames.toArray(new String[0]);
        model.classToIndex = new HashMap<>();
        for (int i = 0; i < classNames.size(); i++) {
            model.classToIndex.put(classNames.get(i), i);
        }
        
        // Converti le strutture degli alberi in DecisionTree
        model.trees = new ArrayList<>();
        for (TreeStructure treeStruct : treeStructures) {
            DecisionTree tree = new DecisionTree(100, 1, model.random);
            tree.setRoot(treeStruct.root);
            model.trees.add(tree);
        }
        
        System.out.println("Modello caricato con successo!");
        return model;
    }
    
    /**
     * Parsing di un singolo albero dal PMML
     */
    private static TreeStructure parseTree(String treeXml, List<String> featureNames) {
        TreeStructure tree = new TreeStructure();
        
        // Trova il primo nodo root (dopo <TreeModel>)
        int firstNodeStart = treeXml.indexOf("<Node");
        if (firstNodeStart == -1) return tree;
        
        tree.root = parseNode(treeXml, firstNodeStart, featureNames);
        return tree;
    }
    
    /**
     * Parsing ricorsivo dei nodi
     */
    private static TreeNode parseNode(String xml, int startPos, List<String> featureNames) {
        // Trova la fine del tag di apertura <Node>
        int nodeOpenEnd = xml.indexOf(">", startPos);
        if (nodeOpenEnd == -1) return null;
        
        // Controlla se è una foglia (ha attributo score)
        String nodeOpenTag = xml.substring(startPos, nodeOpenEnd + 1);
        if (nodeOpenTag.contains("score=\"")) {
            int scoreStart = nodeOpenTag.indexOf("score=\"") + 7;
            int scoreEnd = nodeOpenTag.indexOf("\"", scoreStart);
            String score = nodeOpenTag.substring(scoreStart, scoreEnd);
            return new TreeNode(unescapeXml(score));
        }
        
        // È un nodo interno, cerca i figli
        int currentPos = nodeOpenEnd + 1;
        TreeNode leftChild = null;
        TreeNode rightChild = null;
        int featureIndex = -1;
        double threshold = 0.0;
        
        // Cerca i nodi figli
        int childCount = 0;
        while (childCount < 2) {
            int childNodeStart = xml.indexOf("<Node", currentPos);
            if (childNodeStart == -1) break;
            
            // Trova il SimplePredicate prima di questo nodo
            int predicateStart = xml.lastIndexOf("<SimplePredicate", childNodeStart);
            if (predicateStart > currentPos && predicateStart < childNodeStart) {
                String predicateTag = xml.substring(predicateStart, xml.indexOf(">", predicateStart) + 1);
                
                // Estrai field name
                int fieldStart = predicateTag.indexOf("field=\"") + 7;
                int fieldEnd = predicateTag.indexOf("\"", fieldStart);
                String fieldName = unescapeXml(predicateTag.substring(fieldStart, fieldEnd));
                
                // Estrai operator
                int operatorStart = predicateTag.indexOf("operator=\"") + 10;
                int operatorEnd = predicateTag.indexOf("\"", operatorStart);
                String operator = predicateTag.substring(operatorStart, operatorEnd);
                
                // Estrai value
                int valueStart = predicateTag.indexOf("value=\"") + 7;
                int valueEnd = predicateTag.indexOf("\"", valueStart);
                double value = Double.parseDouble(predicateTag.substring(valueStart, valueEnd));
                
                // Trova l'indice della feature
                featureIndex = featureNames.indexOf(fieldName);
                threshold = value;
                
                // Trova il nodo figlio corrispondente
                int childNodeEnd = findMatchingNodeEnd(xml, childNodeStart);
                TreeNode childNode = parseNode(xml, childNodeStart, featureNames);
                
                if (operator.equals("lessOrEqual")) {
                    leftChild = childNode;
                } else if (operator.equals("greaterThan")) {
                    rightChild = childNode;
                }
                
                currentPos = childNodeEnd;
                childCount++;
            } else {
                break;
            }
        }
        
        if (featureIndex != -1 && leftChild != null && rightChild != null) {
            return new TreeNode(featureIndex, threshold, leftChild, rightChild);
        }
        
        return null;
    }
    
    /**
     * Trova la fine del tag </Node> corrispondente
     */
    private static int findMatchingNodeEnd(String xml, int startPos) {
        int depth = 1;
        int pos = xml.indexOf(">", startPos) + 1;
        
        while (depth > 0 && pos < xml.length()) {
            int nextOpen = xml.indexOf("<Node", pos);
            int nextClose = xml.indexOf("</Node>", pos);
            
            if (nextClose == -1) break;
            
            if (nextOpen != -1 && nextOpen < nextClose) {
                depth++;
                pos = nextOpen + 5;
            } else {
                depth--;
                pos = nextClose + 7;
            }
        }
        
        return pos;
    }
    
    /**
     * Unescape XML entities
     */
    private static String unescapeXml(String text) {
        if (text == null) return "";
        return text.replace("&amp;", "&")
                   .replace("&lt;", "<")
                   .replace("&gt;", ">")
                   .replace("&quot;", "\"")
                   .replace("&apos;", "'");
    }
    
    /**
     * Esporta il modello Random Forest in formato PMML
     */
    public static void exportToPMML(RandomForestClassifier model, String filename) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
            writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            writer.println("<PMML version=\"4.4\" xmlns=\"http://www.dmg.org/PMML-4_4\">");
            writer.println("  <Header copyright=\"Generated by RandomForestFinalModel\">");
            writer.println("    <Application name=\"RandomForestFinalModel\" version=\"1.0\"/>");
            writer.println("    <Timestamp>" + new Date() + "</Timestamp>");
            writer.println("  </Header>");
            
            // Data Dictionary
            writer.println("  <DataDictionary numberOfFields=\"" + (model.getFeatureNames().length + 1) + "\">");
            
            // Features
            for (String featureName : model.getFeatureNames()) {
                writer.println("    <DataField name=\"" + escapeXml(featureName) + "\" optype=\"continuous\" dataType=\"double\"/>");
            }
            
            // Target
            writer.println("    <DataField name=\"Target\" optype=\"categorical\" dataType=\"string\">");
            for (String className : model.getIndexToClass()) {
                writer.println("      <Value value=\"" + escapeXml(className) + "\"/>");
            }
            writer.println("    </DataField>");
            writer.println("  </DataDictionary>");
            
            // Mining Model
            writer.println("  <MiningModel modelName=\"RandomForestModel\" functionName=\"classification\">");
            writer.println("    <MiningSchema>");
            for (String featureName : model.getFeatureNames()) {
                writer.println("      <MiningField name=\"" + escapeXml(featureName) + "\" usageType=\"active\"/>");
            }
            writer.println("      <MiningField name=\"Target\" usageType=\"predicted\"/>");
            writer.println("    </MiningSchema>");
            
            writer.println("    <Segmentation multipleModelMethod=\"majorityVote\">");
            
            // Export each tree
            List<DecisionTree> trees = model.getTrees();
            for (int i = 0; i < trees.size(); i++) {
                writer.println("      <Segment id=\"" + i + "\">");
                writer.println("        <True/>");
                writer.println("        <TreeModel modelName=\"Tree_" + i + "\" functionName=\"classification\" splitCharacteristic=\"binarySplit\">");
                writer.println("          <MiningSchema>");
                for (String featureName : model.getFeatureNames()) {
                    writer.println("            <MiningField name=\"" + escapeXml(featureName) + "\"/>");
                }
                writer.println("            <MiningField name=\"Target\" usageType=\"predicted\"/>");
                writer.println("          </MiningSchema>");
                
                // Export tree structure
                TreeNode root = trees.get(i).getRoot();
                if (root != null) {
                    exportTreeNode(writer, root, model.getFeatureNames(), "          ");
                }
                
                writer.println("        </TreeModel>");
                writer.println("      </Segment>");
            }
            
            writer.println("    </Segmentation>");
            writer.println("  </MiningModel>");
            writer.println("</PMML>");
        }
        
        System.out.println("Modello PMML salvato in: " + filename);
    }
    
    private static void exportTreeNode(PrintWriter writer, TreeNode node, String[] featureNames, String indent) {
        if (node.isLeaf()) {
            writer.println(indent + "<Node score=\"" + escapeXml(node.getPrediction()) + "\">");
            writer.println(indent + "  <True/>");
            writer.println(indent + "</Node>");
        } else {
            writer.println(indent + "<Node>");
            writer.println(indent + "  <True/>");
            
            // Left child (<=)
            if (node.getLeft() != null) {
                writer.println(indent + "  <Node>");
                writer.println(indent + "    <SimplePredicate field=\"" + escapeXml(featureNames[node.getFeatureIndex()]) + 
                             "\" operator=\"lessOrEqual\" value=\"" + node.getThreshold() + "\"/>");
                exportTreeNodeChildren(writer, node.getLeft(), featureNames, indent + "    ");
                writer.println(indent + "  </Node>");
            }
            
            // Right child (>)
            if (node.getRight() != null) {
                writer.println(indent + "  <Node>");
                writer.println(indent + "    <SimplePredicate field=\"" + escapeXml(featureNames[node.getFeatureIndex()]) + 
                             "\" operator=\"greaterThan\" value=\"" + node.getThreshold() + "\"/>");
                exportTreeNodeChildren(writer, node.getRight(), featureNames, indent + "    ");
                writer.println(indent + "  </Node>");
            }
            
            writer.println(indent + "</Node>");
        }
    }
    
    private static void exportTreeNodeChildren(PrintWriter writer, TreeNode node, String[] featureNames, String indent) {
        if (node.isLeaf()) {
            writer.println(indent + "<Node score=\"" + escapeXml(node.getPrediction()) + "\">");
            writer.println(indent + "  <True/>");
            writer.println(indent + "</Node>");
        } else {
            writer.println(indent + "<Node>");
            writer.println(indent + "  <True/>");
            
            if (node.getLeft() != null) {
                writer.println(indent + "  <Node>");
                writer.println(indent + "    <SimplePredicate field=\"" + escapeXml(featureNames[node.getFeatureIndex()]) + 
                             "\" operator=\"lessOrEqual\" value=\"" + node.getThreshold() + "\"/>");
                exportTreeNodeChildren(writer, node.getLeft(), featureNames, indent + "    ");
                writer.println(indent + "  </Node>");
            }
            
            if (node.getRight() != null) {
                writer.println(indent + "  <Node>");
                writer.println(indent + "    <SimplePredicate field=\"" + escapeXml(featureNames[node.getFeatureIndex()]) + 
                             "\" operator=\"greaterThan\" value=\"" + node.getThreshold() + "\"/>");
                exportTreeNodeChildren(writer, node.getRight(), featureNames, indent + "    ");
                writer.println(indent + "  </Node>");
            }
            
            writer.println(indent + "</Node>");
        }
    }
    
    private static String escapeXml(String text) {
        if (text == null) return "";
        return text.replace("&", "&amp;")
                   .replace("<", "&lt;")
                   .replace(">", "&gt;")
                   .replace("\"", "&quot;")
                   .replace("'", "&apos;");
    }
}

/* -------------------------
   PredictionResult, Instance, DataSet, EvaluationStats
   ------------------------- */
class PredictionResult {
    private String prediction;
    private double confidence;

    public PredictionResult(String prediction, double confidence) {
        this.prediction = prediction;
        this.confidence = confidence;
    }

    public String getPrediction() { return prediction; }
    public double getConfidence() { return confidence; }
}

class Instance {
    private String ds;
    private String id;
    private String smiles;
    private String target;
    private double[] features;

    public Instance(String ds, String id, String smiles, String target, double[] features) {
        this.ds = ds;
        this.id = id;
        this.smiles = smiles;
        this.target = target;
        this.features = features;
    }

    public String getDs() { return ds; }
    public String getId() { return id; }
    public String getSmiles() { return smiles; }
    public String getTarget() { return target; }
    public double[] getFeatures() { return features; }
    public double getFeatureValue(int index) { return features[index]; }
}

class DataSet {
    private List<Instance> instances;
    private String[] featureNames;
    private Map<String, Integer> classMapping;

    public DataSet(List<Instance> instances, String[] featureNames, Map<String, Integer> classMapping) {
        this.instances = instances;
        this.featureNames = featureNames;
        this.classMapping = classMapping;
    }

    public List<Instance> getInstances() { return instances; }
    public String[] getFeatureNames() { return featureNames; }
    public Map<String, Integer> getClassMapping() { return classMapping; }
}

class EvaluationStats {
    private int tp, fp, tn, fn, outOfDomain;
    private double balancedAccuracy, accuracy, sensitivity, specificity, ppv, npv, mcc;

    public EvaluationStats(List<String> actual, List<String> predicted) {
        this.outOfDomain = 0;
        this.tp = this.fp = this.tn = this.fn = 0;

        Set<String> labels = new HashSet<>();
        labels.addAll(actual);
        labels.addAll(predicted);
        List<String> labs = new ArrayList<>(labels);
        if (labs.size() < 2) {
            this.balancedAccuracy = this.accuracy = this.sensitivity = this.specificity =
                this.ppv = this.npv = this.mcc = 1.0;
            return;
        }
        String pos = labs.get(0);
        String neg = labs.get(1);

        for (int i = 0; i < actual.size(); i++) {
            String a = actual.get(i);
            String p = predicted.get(i);
            if (a.equals(pos) && p.equals(pos)) tp++;
            if (a.equals(neg) && p.equals(pos)) fp++;
            if (a.equals(neg) && p.equals(neg)) tn++;
            if (a.equals(pos) && p.equals(neg)) fn++;
        }

        double sens = tp + fn == 0 ? 0.0 : (double) tp / (tp + fn);
        double spec = tn + fp == 0 ? 0.0 : (double) tn / (tn + fp);
        double acc = (double) (tp + tn) / (tp + tn + fp + fn + 1e-12);
        double ppv_ = tp + fp == 0 ? 0.0 : (double) tp / (tp + fp);
        double npv_ = tn + fn == 0 ? 0.0 : (double) tn / (tn + fn);

        this.sensitivity = sens;
        this.specificity = spec;
        this.accuracy = acc;
        this.balancedAccuracy = (sens + spec) / 2.0;
        this.ppv = ppv_;
        this.npv = npv_;
        double denom = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
        this.mcc = denom == 0 ? 0.0 : ((double) tp * tn - fp * fn) / denom;
    }

    public int getTp() { return tp; }
    public int getFp() { return fp; }
    public int getTn() { return tn; }
    public int getFn() { return fn; }
    public int getOutOfDomain() { return outOfDomain; }
    public double getBalancedAccuracy() { return balancedAccuracy; }
    public double getAccuracy() { return accuracy; }
    public double getSensitivity() { return sensitivity; }
    public double getSpecificity() { return specificity; }
    public double getPpv() { return ppv; }
    public double getNpv() { return npv; }
    public double getMcc() { return mcc; }
}

/* -------------------------
   Programma principale
   ------------------------- */
public class RandomForestFinalModel {

    private static final String[] FEATURE_NAMES = {
        "MAXDP", "IDE", "SIC1", "SpMaxA_D/Dt", "SpMax_B(s)",
        "GATS7i", "SpMax2_Bh(s)", "P_VSA_ppp_ar", "P_VSA_ppp_con",
        "Eta_beta_A", "SM05_EA(bo)", "SM02_EA(dm)", "Eig12_AEA(bo)", "SaaCH"
    };

    private static final int OPTIMAL_TREES = 500;
    private static final int OPTIMAL_MAX_DEPTH = 10;
    private static final int OPTIMAL_MIN_NODE_SIZE = 2;

    public static void main(String[] args) {
        try {
            System.out.println("=== Random Forest Final Model ===");
            System.out.println("Usando i parametri ottimali trovati:");
            System.out.println("- Trees: " + OPTIMAL_TREES);
            System.out.println("- Max Depth: " + OPTIMAL_MAX_DEPTH);
            System.out.println("- Min Node Size: " + OPTIMAL_MIN_NODE_SIZE);
            System.out.println();

            System.out.println("1. Caricamento dati...");
            DataSet trainingSet = loadDataSet("TS.csv");
            DataSet calibrationSet = loadDataSet("CS.csv");
            DataSet validationSet = loadDataSet("VS.csv");

            System.out.println("\n2. Training del modello...");
            RandomForestClassifier model = new RandomForestClassifier(
                OPTIMAL_TREES, OPTIMAL_MAX_DEPTH, OPTIMAL_MIN_NODE_SIZE);
            model.train(trainingSet);

            // Salva il modello in formato .ser (serializzazione Java)
            try {
                saveModel(model, "RandomForestModel.ser");
                System.out.println("Modello serializzato salvato in: RandomForestModel.ser");
            } catch (IOException e) {
                System.err.println("Attenzione: impossibile salvare il modello .ser: " + e.getMessage());
            }

            // Salva il modello in formato PMML (XML standard)
            try {
                PMMLExporter.exportToPMML(model, "RandomForestModel.pmml");
                System.out.println("Modello PMML salvato in: RandomForestModel.pmml");
            } catch (IOException e) {
                System.err.println("Attenzione: impossibile salvare il modello PMML: " + e.getMessage());
            }

            // TEST: Ricarica il modello dal PMML per verificare la correttezza
            System.out.println("\n=== TEST: Ricaricamento modello da PMML ===");
            try {
                RandomForestClassifier loadedModel = PMMLExporter.importFromPMML("RandomForestModel.pmml");
                
                // Verifica predizioni su un subset del training set
                System.out.println("Verifica predizioni modello caricato...");
                List<Instance> testInstances = trainingSet.getInstances().subList(0, Math.min(10, trainingSet.getInstances().size()));
                
                System.out.println("Confronto predizioni (prime 10 istanze):");
                System.out.println("ID | Originale | Caricato | Match");
                System.out.println("----------------------------------------");
                
                boolean allMatch = true;
                for (Instance inst : testInstances) {
                    PredictionResult origPred = model.predict(inst);
                    PredictionResult loadedPred = loadedModel.predict(inst);
                    boolean match = origPred.getPrediction().equals(loadedPred.getPrediction());
                    allMatch = allMatch && match;
                    
                    System.out.printf("%s | %s | %s | %s%n",
                        inst.getId(),
                        origPred.getPrediction(),
                        loadedPred.getPrediction(),
                        match ? "✓" : "✗");
                }
                
                if (allMatch) {
                    System.out.println("\n✓ Modello PMML caricato correttamente! Tutte le predizioni coincidono.");
                } else {
                    System.out.println("\n✗ Attenzione: alcune predizioni differiscono.");
                }
                
            } catch (Exception e) {
                System.err.println("Errore nel test di ricaricamento: " + e.getMessage());
                e.printStackTrace();
            }

            System.out.println("\n3. Predizioni...");
            System.out.println("Predizione Training Set:");
            List<PredictionResult> tsPredictions = model.predict(trainingSet);

            System.out.println("Predizione Calibration Set:");
            List<PredictionResult> csPredictions = model.predict(calibrationSet);

            System.out.println("Predizione Validation Set:");
            List<PredictionResult> vsPredictions = model.predict(validationSet);

            System.out.println("\n4. Salvataggio predizioni...");
            savePredictions(trainingSet, tsPredictions, calibrationSet, csPredictions,
                          validationSet, vsPredictions);

            System.out.println("5. Calcolo statistiche...");
            saveStatistics(trainingSet, tsPredictions, calibrationSet, csPredictions,
                          validationSet, vsPredictions);

            System.out.println("\n=== Riassunto Performance ===");
            printPerformanceSummary(trainingSet, tsPredictions, "Training Set");
            printPerformanceSummary(calibrationSet, csPredictions, "Calibration Set");
            printPerformanceSummary(validationSet, vsPredictions, "Validation Set");

        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("Errore durante l'esecuzione: " + e.getMessage());
        }
    }

    private static void saveModel(RandomForestClassifier model, String filename) throws IOException {
        try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filename))) {
            oos.writeObject(model);
        }
    }

    private static RandomForestClassifier loadModel(String filename) throws IOException, ClassNotFoundException {
        try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(filename))) {
            return (RandomForestClassifier) ois.readObject();
        }
    }

    private static DataSet loadDataSet(String filename) throws IOException {
        List<Instance> instances = new ArrayList<>();
        Map<String, Integer> classMapping = new HashMap<>();

        System.out.println("Caricamento file: " + filename);

        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IOException("File vuoto: " + filename);
            }

            String[] headerParts = headerLine.split(";");
            Map<String, Integer> headerIndex = new HashMap<>();
            for (int i = 0; i < headerParts.length; i++) {
                headerIndex.put(headerParts[i].trim(), i);
            }

            int idxDS = headerIndex.getOrDefault("DS", 0);
            int idxID = headerIndex.getOrDefault("ID", 1);
            int idxSMILES = headerIndex.getOrDefault("SMILES", 2);
            int idxTarget = headerIndex.getOrDefault("Target", 3);

            int[] featureIndices = new int[FEATURE_NAMES.length];
            boolean allFound = true;
            for (int i = 0; i < FEATURE_NAMES.length; i++) {
                String fname = FEATURE_NAMES[i];
                if (headerIndex.containsKey(fname)) {
                    featureIndices[i] = headerIndex.get(fname);
                } else {
                    allFound = false;
                    break;
                }
            }

            int lineNumber = 1;
            String line;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                String[] parts = line.split(";", -1);

                if (!allFound) {
                    if (parts.length < 4 + FEATURE_NAMES.length) {
                        System.out.println("Riga " + lineNumber + " ignorata - troppe poche colonne: " + parts.length);
                        continue;
                    }
                }

                String ds = parts.length > idxDS ? parts[idxDS].trim() : "";
                String id = parts.length > idxID ? parts[idxID].trim() : "";
                String smiles = parts.length > idxSMILES ? parts[idxSMILES].trim() : "";
                String target = parts.length > idxTarget ? parts[idxTarget].trim() : "";

                if (!classMapping.containsKey(target)) {
                    classMapping.put(target, classMapping.size());
                }

                double[] features = new double[FEATURE_NAMES.length];

                if (allFound) {
                    for (int i = 0; i < FEATURE_NAMES.length; i++) {
                        int col = featureIndices[i];
                        try {
                            String featureValue = col < parts.length ? parts[col].trim().replace(",", ".") : "";
                            features[i] = featureValue.isEmpty() ? 0.0 : Double.parseDouble(featureValue);
                        } catch (NumberFormatException e) {
                            features[i] = 0.0;
                        }
                    }
                } else {
                    for (int i = 0; i < FEATURE_NAMES.length; i++) {
                        try {
                            String featureValue = parts[4 + i].trim().replace(",", ".");
                            features[i] = featureValue.isEmpty() ? 0.0 : Double.parseDouble(featureValue);
                        } catch (Exception e) {
                            features[i] = 0.0;
                        }
                    }
                }

                instances.add(new Instance(ds, id, smiles, target, features));
            }
        }

        Map<String, Long> classDistribution = instances.stream()
            .collect(Collectors.groupingBy(Instance::getTarget, Collectors.counting()));
        System.out.println("- Istanze caricate: " + instances.size());
        System.out.println("- Distribuzione classi: " + classDistribution);

        return new DataSet(instances, FEATURE_NAMES, classMapping);
    }

    private static void savePredictions(DataSet trainingSet, List<PredictionResult> tsPredictions,
                                      DataSet calibrationSet, List<PredictionResult> csPredictions,
                                      DataSet validationSet, List<PredictionResult> vsPredictions) throws IOException {

        try (PrintWriter writer = new PrintWriter(new FileWriter("predictions.csv"))) {
            writer.println("DS;ID;SMILES;Target;Prediction;Confidence");

            saveDatasetPredictions(writer, trainingSet, tsPredictions);
            saveDatasetPredictions(writer, calibrationSet, csPredictions);
            saveDatasetPredictions(writer, validationSet, vsPredictions);
        }

        System.out.println("Predizioni salvate in: predictions.csv");
    }

    private static void saveDatasetPredictions(PrintWriter writer, DataSet dataSet,
                                             List<PredictionResult> predictions) {
        List<Instance> instances = dataSet.getInstances();
        for (int i = 0; i < instances.size(); i++) {
            Instance instance = instances.get(i);
            PredictionResult prediction = predictions.get(i);

            writer.printf("%s;%s;%s;%s;%s;%.4f%n",
                instance.getDs(), instance.getId(), instance.getSmiles(),
                instance.getTarget(), prediction.getPrediction(),
                prediction.getConfidence());
        }
    }

    private static void saveStatistics(DataSet trainingSet, List<PredictionResult> tsPredictions,
                                     DataSet calibrationSet, List<PredictionResult> csPredictions,
                                     DataSet validationSet, List<PredictionResult> vsPredictions) throws IOException {

        try (PrintWriter writer = new PrintWriter(new FileWriter("stats.csv"))) {
            writer.println("Dataset;Balanced_Accuracy;Accuracy;Sensitivity;Specificity;PPV;NPV;MCC;TP;FP;TN;FN;Out_of_Domain");

            EvaluationStats tsStats = calculateStats(trainingSet, tsPredictions);
            writeStats(writer, "TS", tsStats);

            EvaluationStats csStats = calculateStats(calibrationSet, csPredictions);
            writeStats(writer, "CS", csStats);

            EvaluationStats vsStats = calculateStats(validationSet, vsPredictions);
            writeStats(writer, "VS", vsStats);
        }
    }

    private static EvaluationStats calculateStats(DataSet dataSet, List<PredictionResult> predictions) {
        List<String> actual = dataSet.getInstances().stream()
            .map(Instance::getTarget)
            .collect(Collectors.toList());

        List<String> predicted = predictions.stream()
            .map(PredictionResult::getPrediction)
            .collect(Collectors.toList());

        return new EvaluationStats(actual, predicted);
    }

    private static void writeStats(PrintWriter writer, String datasetName, EvaluationStats stats) {
        writer.printf("%s;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%d;%d;%d;%d;%d%n",
            datasetName, stats.getBalancedAccuracy(), stats.getAccuracy(),
            stats.getSensitivity(), stats.getSpecificity(), stats.getPpv(),
            stats.getNpv(), stats.getMcc(), stats.getTp(), stats.getFp(),
            stats.getTn(), stats.getFn(), stats.getOutOfDomain());
    }

    private static void printPerformanceSummary(DataSet dataSet, List<PredictionResult> predictions, String setName) {
        EvaluationStats stats = calculateStats(dataSet, predictions);

        System.out.println(setName + ":");
        System.out.printf("  Accuracy: %.4f, Balanced Accuracy: %.4f%n",
            stats.getAccuracy(), stats.getBalancedAccuracy());
        System.out.printf("  Sensitivity: %.4f, Specificity: %.4f%n",
            stats.getSensitivity(), stats.getSpecificity());
        System.out.printf("  PPV: %.4f, NPV: %.4f%n",
            stats.getPpv(), stats.getNpv());
        System.out.printf("  MCC: %.4f%n", stats.getMcc());
        System.out.printf("  TP: %d, FP: %d, TN: %d, FN: %d%n",
            stats.getTp(), stats.getFp(), stats.getTn(), stats.getFn());
    }
}