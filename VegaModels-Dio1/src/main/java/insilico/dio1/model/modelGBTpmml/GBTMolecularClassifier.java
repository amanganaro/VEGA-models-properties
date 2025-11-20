package insilico.dio1.model.modelGBTpmml;

import insilico.core.exception.InitFailureException;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// Modello GBT ottimizzato con parametri fissi + export/import PMML
public class GBTMolecularClassifier {

    private final static String PMML_PATH = "/data/dio1_gbt_model.pmml";
    private final GradientBoostedTrees gbt;
    private final DioDescriptors descriptors;

    public GBTMolecularClassifier() throws InitFailureException {
        try {
            this.gbt = new GradientBoostedTrees(PMML_PATH);
            this.descriptors = new DioDescriptors();
        } catch (Throwable e) {
            throw new InitFailureException(e);
        }
    }

    public GBTResult predict(String SMILES)  {
        GBTResult res = new GBTResult();

        try {
            double[] desc = descriptors.Scale(descriptors.Calculate(SMILES));

            // round - original model trained on values with 4 significant digits
            for (int i=0; i<desc.length; i++)
                desc[i] = Math.round( desc[i] * 10000.0 ) / 10000.0;

            res.probability = this.gbt.predictProbability(desc);
            res.prediction = this.gbt.predictClass(desc);
            res.confidence = Math.abs(res.probability - 0.5) * 2.0;
        } catch (Exception e) {
            // do nothing - res already initialized to missing value
        }

        return res;
    }

    // Struttura dati per rappresentare un campione
    public class Sample {
        String ds, id, smiles;
        int target;
        double[] features;
        int prediction;
        double probability;
        double confidence;
        
        public Sample(String ds, String id, String smiles, int target, double[] features) {
            this.ds = ds;
            this.id = id;
            this.smiles = smiles;
            this.target = target;
            this.features = features;
        }
    }
    
    // Classe per rappresentare un nodo dell'albero
    public class TreeNode {
        int featureIndex = -1;
        double threshold;
        double value;
        TreeNode left, right;
        boolean isLeaf = false;
        String nodeId; // Per PMML
        
        public TreeNode(double value) {
            this.value = value;
            this.isLeaf = true;
        }
        
        public TreeNode(int featureIndex, double threshold) {
            this.featureIndex = featureIndex;
            this.threshold = threshold;
        }
    }
    
    // Implementazione del Gradient Boosted Trees con parametri fissi
    public class GradientBoostedTrees {
        List<TreeNode> trees;
        int maxDepth = 5;
        int minNodeSize = 3;
        int numTrees = 500;
        double learningRate = 0.1;
        double initialPrediction = 0.0;
        
        public GradientBoostedTrees() {
            this.trees = new ArrayList<>();
        }
        
        public void train(List<Sample> samples) {
            System.out.println("Training GBT con parametri: trees=500, maxDepth=5, minNodeSize=3");
            
            initialPrediction = samples.stream().mapToInt(s -> s.target).average().orElse(0.0);
            double[] predictions = new double[samples.size()];
            Arrays.fill(predictions, initialPrediction);
            
            for (int t = 0; t < numTrees; t++) {
                if (t % 100 == 0) {
                    System.out.printf("Training albero %d/%d%n", t + 1, numTrees);
                }
                
                double[] residuals = new double[samples.size()];
                for (int i = 0; i < samples.size(); i++) {
                    double probability = sigmoid(predictions[i]);
                    residuals[i] = samples.get(i).target - probability;
                }
                
                TreeNode tree = buildTree(samples, residuals, 0);
                trees.add(tree);
                
                for (int i = 0; i < samples.size(); i++) {
                    double treePredict = predictTree(tree, samples.get(i).features);
                    predictions[i] += learningRate * treePredict;
                }
            }
            System.out.println("Training completato!");
        }
        
        private TreeNode buildTree(List<Sample> samples, double[] residuals, int depth) {
            if (depth >= maxDepth || samples.size() < minNodeSize) {
                double avgResidual = 0;
                for (int i = 0; i < residuals.length; i++) {
                    avgResidual += residuals[i];
                }
                return new TreeNode(avgResidual / residuals.length);
            }
            
            int bestFeature = -1;
            double bestThreshold = 0;
            double bestGain = Double.NEGATIVE_INFINITY;
            int numFeatures = samples.get(0).features.length;
            
            for (int f = 0; f < numFeatures; f++) {
                Set<Double> uniqueValues = new HashSet<>();
                for (Sample s : samples) {
                    uniqueValues.add(s.features[f]);
                }
                
                for (double threshold : uniqueValues) {
                    double gain = calculateGain(samples, residuals, f, threshold);
                    if (gain > bestGain) {
                        bestGain = gain;
                        bestFeature = f;
                        bestThreshold = threshold;
                    }
                }
            }
            
            if (bestFeature == -1) {
                double avgResidual = 0;
                for (double r : residuals) avgResidual += r;
                return new TreeNode(avgResidual / residuals.length);
            }
            
            List<Sample> leftSamples = new ArrayList<>();
            List<Sample> rightSamples = new ArrayList<>();
            List<Double> leftResiduals = new ArrayList<>();
            List<Double> rightResiduals = new ArrayList<>();
            
            for (int i = 0; i < samples.size(); i++) {
                if (samples.get(i).features[bestFeature] <= bestThreshold) {
                    leftSamples.add(samples.get(i));
                    leftResiduals.add(residuals[i]);
                } else {
                    rightSamples.add(samples.get(i));
                    rightResiduals.add(residuals[i]);
                }
            }
            
            TreeNode node = new TreeNode(bestFeature, bestThreshold);
            if (!leftSamples.isEmpty()) {
                node.left = buildTree(leftSamples, leftResiduals.stream().mapToDouble(Double::doubleValue).toArray(), depth + 1);
            }
            if (!rightSamples.isEmpty()) {
                node.right = buildTree(rightSamples, rightResiduals.stream().mapToDouble(Double::doubleValue).toArray(), depth + 1);
            }
            
            return node;
        }
        
        private double calculateGain(List<Sample> samples, double[] residuals, int feature, double threshold) {
            double totalVariance = calculateVariance(residuals);
            
            List<Double> leftResiduals = new ArrayList<>();
            List<Double> rightResiduals = new ArrayList<>();
            
            for (int i = 0; i < samples.size(); i++) {
                if (samples.get(i).features[feature] <= threshold) {
                    leftResiduals.add(residuals[i]);
                } else {
                    rightResiduals.add(residuals[i]);
                }
            }
            
            if (leftResiduals.isEmpty() || rightResiduals.isEmpty()) {
                return Double.NEGATIVE_INFINITY;
            }
            
            double leftVar = calculateVariance(leftResiduals.stream().mapToDouble(Double::doubleValue).toArray());
            double rightVar = calculateVariance(rightResiduals.stream().mapToDouble(Double::doubleValue).toArray());
            double weightedVar = (leftResiduals.size() * leftVar + rightResiduals.size() * rightVar) / residuals.length;
            
            return totalVariance - weightedVar;
        }
        
        private double calculateVariance(double[] values) {
            double mean = Arrays.stream(values).average().orElse(0.0);
            return Arrays.stream(values).map(x -> Math.pow(x - mean, 2)).average().orElse(0.0);
        }
        
        public double predictProbability(double[] features) {
            double prediction = initialPrediction;
            for (TreeNode tree : trees) {
                prediction += learningRate * predictTree(tree, features);
            }
            return sigmoid(prediction);
        }
        
        public int predictClass(double[] features) {
            double probability = predictProbability(features);
            return probability >= 0.5 ? 1 : 0;
        }
        
        private double predictTree(TreeNode node, double[] features) {
            if (node.isLeaf) {
                return node.value;
            }
            
            if (features[node.featureIndex] <= node.threshold) {
                return node.left != null ? predictTree(node.left, features) : 0;
            } else {
                return node.right != null ? predictTree(node.right, features) : 0;
            }
        }
        
        private double sigmoid(double x) {
            return 1.0 / (1.0 + Math.exp(-x));
        }
        
        // ===== EXPORT PMML =====
        public void exportToPMML(String filename) throws IOException {
            System.out.println("Esportazione modello in PMML...");
            PrintWriter writer = new PrintWriter(new FileWriter(filename));
            
            writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            writer.println("<PMML version=\"4.4\" xmlns=\"http://www.dmg.org/PMML-4_4\">");
            writer.println("  <Header>");
            writer.println("    <Application name=\"GBTMolecularClassifier\" version=\"1.0\"/>");
            writer.println("  </Header>");
            writer.println("  <DataDictionary numberOfFields=\"18\">");
            
            for (int i = 0; i < 17; i++) {
                writer.println("    <DataField name=\"feature" + i + "\" optype=\"continuous\" dataType=\"double\"/>");
            }
            writer.println("    <DataField name=\"target\" optype=\"categorical\" dataType=\"integer\">");
            writer.println("      <Value value=\"0\"/>");
            writer.println("      <Value value=\"1\"/>");
            writer.println("    </DataField>");
            writer.println("  </DataDictionary>");
            
            writer.println("  <MiningModel functionName=\"classification\">");
            writer.println("    <MiningSchema>");
            for (int i = 0; i < 17; i++) {
                writer.println("      <MiningField name=\"feature" + i + "\" usageType=\"active\"/>");
            }
            writer.println("      <MiningField name=\"target\" usageType=\"target\"/>");
            writer.println("    </MiningSchema>");
            
            writer.println("    <Segmentation multipleModelMethod=\"weightedSum\">");
            writer.println("      <Segment id=\"initialValue\">");
            writer.println("        <True/>");
            writer.println("        <Regression>");
            writer.println("          <RegressionTable intercept=\"" + initialPrediction + "\"/>");
            writer.println("        </Regression>");
            writer.println("      </Segment>");
            
            for (int i = 0; i < trees.size(); i++) {
                writer.println("      <Segment id=\"tree" + i + "\" weight=\"" + learningRate + "\">");
                writer.println("        <True/>");
                writer.println("        <TreeModel functionName=\"regression\" splitCharacteristic=\"binarySplit\">");
                writer.println("          <MiningSchema>");
                for (int f = 0; f < 17; f++) {
                    writer.println("            <MiningField name=\"feature" + f + "\"/>");
                }
                writer.println("          </MiningSchema>");
                
                int[] nodeCounter = {0};
                exportTreeNode(writer, trees.get(i), nodeCounter, "          ");
                
                writer.println("        </TreeModel>");
                writer.println("      </Segment>");
            }
            
            writer.println("    </Segmentation>");
            writer.println("  </MiningModel>");
            writer.println("</PMML>");
            
            writer.close();
            System.out.println("Modello esportato in: " + filename);
        }
        
        private void exportTreeNode(PrintWriter writer, TreeNode node, int[] counter, String indent) {
            String nodeId = "n" + counter[0]++;
            node.nodeId = nodeId;
            
            if (node.isLeaf) {
                writer.println(indent + "<Node id=\"" + nodeId + "\" score=\"" + node.value + "\"/>");
            } else {
                writer.println(indent + "<Node id=\"" + nodeId + "\">");
                writer.println(indent + "  <SimplePredicate field=\"feature" + node.featureIndex + 
                             "\" operator=\"lessOrEqual\" value=\"" + node.threshold + "\"/>");
                
                if (node.left != null) {
                    exportTreeNode(writer, node.left, counter, indent + "  ");
                }
                if (node.right != null) {
                    exportTreeNode(writer, node.right, counter, indent + "  ");
                }
                
                writer.println(indent + "</Node>");
            }
        }
        
        // ===== IMPORT PMML =====
        // trasformato in costruttore della classe, overload
        public GradientBoostedTrees(String filename) throws IOException {
            super();
            this.trees = new ArrayList<>();

            URL DataSource = GBTMolecularClassifier.class.getResource(filename);
            BufferedReader reader = new BufferedReader(new InputStreamReader(DataSource.openStream()));

            StringBuilder content = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                content.append(line).append("\n");
            }
            reader.close();
            
            String pmml = content.toString();
            
            // Estrai initialPrediction
            Pattern initialPattern = Pattern.compile("<Segment id=\"initialValue\">.*?<RegressionTable intercept=\"([^\"]+)\"", 
                                                    Pattern.DOTALL);
            Matcher initialMatcher = initialPattern.matcher(pmml);
            if (initialMatcher.find()) {
                this.initialPrediction = Double.parseDouble(initialMatcher.group(1));
            }
            
            // Estrai learning rate dal primo segment con weight
            Pattern weightPattern = Pattern.compile("<Segment id=\"tree\\d+\" weight=\"([^\"]+)\"");
            Matcher weightMatcher = weightPattern.matcher(pmml);
            if (weightMatcher.find()) {
                this.learningRate = Double.parseDouble(weightMatcher.group(1));
            }
            
            // Estrai ogni albero
            Pattern treePattern = Pattern.compile("<Segment id=\"tree\\d+\"[^>]*>.*?</Segment>", Pattern.DOTALL);
            Matcher treeMatcher = treePattern.matcher(pmml);
            
            while (treeMatcher.find()) {
                String treeXml = treeMatcher.group();
                TreeNode tree = parseTreeFromXML(treeXml);
                this.trees.add(tree);
            }

            this.numTrees = this.trees.size();
        }
        
        private TreeNode parseTreeFromXML(String xml) {

//            String treeXml = getClosedTag(xml,"node", 0).xml;
            String treeXml = XmlTagUtils.getClosedTag(xml, "node").get("result");
            if (treeXml.isEmpty())
                return new TreeNode(0.0);
            return parseNode(treeXml);

            // Trova il primo nodo root
//            Pattern rootPattern = Pattern.compile("<Node id=\"([^\"]+)\">(.*?)</Node>", Pattern.DOTALL);
//            Matcher rootMatcher = rootPattern.matcher(xml);
//
//            if (rootMatcher.find()) {
//                return parseNode(rootMatcher.group(0));
//            }
//
//            return new TreeNode(0.0);
        }

//        private class TagResult {
//            String xml = "";
//            int finalIdx = 0;
//        }
//
//        // recupera tutto il contenuto del tag
//        private TagResult getClosedTag(String xml, String tagName, int startIdx) {
//            String res = "";
//
//            int idx = xml.toLowerCase().indexOf("<" + tagName.toLowerCase(), startIdx);
//            if (idx == -1)
//                return new TagResult();
////            for (int i=(idx+1+tagName.length()) ; i<xml.length(); i++) {
////                if (xml.charAt(i)=='>') {
////                    idx = i;
////                    break;
////                }
////            }
//
//            int nOpen = 0;
//            for (int i=idx ; i<xml.length(); i++) {
//
//                if (xml.charAt(i) == '<') {
//
//                    if (xml.toLowerCase().indexOf("<" + tagName.toLowerCase(), i) == i) {
//                        // aperto un tag
//
//                        // controlla se si chiude senza tag chiusura
//                        int j = i + 1;
//                        boolean closed = false;
//                        while (j < xml.length()) {
//                            if (xml.charAt(j) == '>')
//                                break;
//                            if ((xml.charAt(j) == '/') && (xml.charAt(j + 1) == '>')) {
//                                closed = true;
//                                break;
//                            }
//                            j++;
//                        }
//
//                        if (!closed)
//                            nOpen++;
//
//                        res += xml.charAt(i);
//
//                    } else if (xml.toLowerCase().indexOf("</" + tagName.toLowerCase(), i) == i) {
//
//                        if (nOpen == 1) {
//                            TagResult r = new TagResult();
//                            String closing = "</" + tagName + ">";
//                            r.xml =  res + closing;
//                            r.finalIdx = i + closing.length();
//                            return r;
//                        } else {
//                            res += xml.charAt(i);
//                            nOpen--;
//                        }
//
//                    } else {
//                        res += xml.charAt(i);
//                    }
//
//                } else {
//                    res += xml.charAt(i);
//                }
//            }
//
//            return new TagResult();
//        }

        private String getFirstTag(String xml, int startIdx) {
            String res = "";
            boolean open = false;
            for (int i=startIdx; i<xml.length(); i++) {
                if (!open) {
                    if (xml.charAt(i) == '<') {
                        open = true;
                        continue;
                    }
                }
                else {
                    if (xml.charAt(i) == '>')
                        break;
                    else
                        res += xml.charAt(i);
                }
            }
            return res;
        }

        private  TreeNode parseNode(String nodeXml) {

            ArrayList<String> tags = XmlTagUtils.divideTags(nodeXml);
            if (tags.get(0).toLowerCase().contains("score=")) {

            // Check se è foglia (ha score)
            Pattern scorePattern = Pattern.compile("<Node id=\"[^\"]+\" score=\"([^\"]+)\"/>");
            Matcher scoreMatcher = scorePattern.matcher(nodeXml);

            if (scoreMatcher.find()) {
                    return new TreeNode(Double.parseDouble(scoreMatcher.group(1)));
                }
            }
            
            // Altrimenti è nodo interno - estrai predicate
            Pattern predicatePattern = Pattern.compile("<SimplePredicate field=\"feature(\\d+)\" operator=\"lessOrEqual\" value=\"([^\"]+)\"/>");
            Matcher predicateMatcher = predicatePattern.matcher(nodeXml);
            
            if (!predicateMatcher.find()) {
                return new TreeNode(0.0);
            }
            
            int featureIndex = Integer.parseInt(predicateMatcher.group(1));
            double threshold = Double.parseDouble(predicateMatcher.group(2));
            
            TreeNode node = new TreeNode(featureIndex, threshold);

            HashMap<String, String> res = XmlTagUtils.getClosedTag(nodeXml.substring(2), "node");
            String childLeft = res.get("result");
            if (!childLeft.isEmpty()) {
                node.left = parseNode(childLeft);
                String childRight = XmlTagUtils.getClosedTag(res.get("remaining"), "node").get("result");
                if (!childRight.isEmpty())
                    node.right = parseNode(childRight);
            }

//            int idx = getFirstTag(nodeXml, 0).length();
//            TagResult left = getClosedTag(nodeXml, "node", idx);
//            String childLeft = left.xml;
//            if (!childLeft.isEmpty()) {
//                node.left = parseNode(childLeft);
//                String childRight = getClosedTag(nodeXml, "node", left.finalIdx).xml;
//                if (!childRight.isEmpty())
//                    node.right = parseNode(childRight);
//            }



//            // Trova nodi figli
//            Pattern childPattern = Pattern.compile("<Node id=\"[^\"]+\"[^>]*>.*?</Node>", Pattern.DOTALL);
//            Matcher childMatcher = childPattern.matcher(nodeXml);
//
//            List<String> children = new ArrayList<>();
//            int lastEnd = 0;
//            while (childMatcher.find()) {
//                if (childMatcher.start() > lastEnd) {
//                    children.add(childMatcher.group());
//                    lastEnd = childMatcher.end();
//                }
//            }
//
//            if (children.size() >= 1) {
//                node.left = parseNode(children.get(0));
//            }
//            if (children.size() >= 2) {
//                node.right = parseNode(children.get(1));
//            }
            
            return node;
        }
    }
    
    public class ClassificationStats {
        int tp, fp, tn, fn;
        double accuracy, balancedAccuracy, sensitivity, specificity, ppv, npv, mcc;
        
        public ClassificationStats(List<Sample> samples) {
            calculate(samples);
        }
        
        private void calculate(List<Sample> samples) {
            tp = fp = tn = fn = 0;
            
            for (Sample s : samples) {
                if (s.prediction == 1 && s.target == 1) tp++;
                else if (s.prediction == 1 && s.target == 0) fp++;
                else if (s.prediction == 0 && s.target == 0) tn++;
                else if (s.prediction == 0 && s.target == 1) fn++;
            }
            
            accuracy = (double)(tp + tn) / (tp + fp + tn + fn);
            sensitivity = tp > 0 ? (double)tp / (tp + fn) : 0;
            specificity = tn > 0 ? (double)tn / (tn + fp) : 0;
            balancedAccuracy = (sensitivity + specificity) / 2.0;
            ppv = tp > 0 ? (double)tp / (tp + fp) : 0;
            npv = tn > 0 ? (double)tn / (tn + fn) : 0;
            
            double denominator = Math.sqrt((long)(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
            mcc = denominator > 0 ? (double)(tp * tn - fp * fn) / denominator : 0;
        }
    }
    
    public List<Sample> readCSV(String filename) throws IOException {
        List<Sample> samples = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String headerLine = br.readLine();
        
        System.out.println("\n=== Lettura file: " + filename + " ===");
        System.out.println("Header: " + headerLine);
        
        String[] headerParts = headerLine.split(",");
        System.out.println("Numero colonne header: " + headerParts.length);
        
        int lineNumber = 1;
        int skippedLines = 0;
        int loadedLines = 0;
        
        String line;
        while ((line = br.readLine()) != null) {
            lineNumber++;
            line = line.trim();
            
            if (line.isEmpty()) {
                continue;
            }
            
            String[] parts = line.split(",");
            
            if (parts.length != 21) {
                System.err.println("\n=== ERRORE Riga " + lineNumber + " ===");
                System.err.println("Colonne trovate: " + parts.length + " (attese: 21)");
                System.err.println("Riga completa: " + line);
                System.err.println("Prime colonne:");
                for (int i = 0; i < Math.min(10, parts.length); i++) {
                    System.err.println("  [" + i + "]: '" + parts[i] + "'");
                }
                skippedLines++;
                continue;
            }
            
            try {
                String ds = parts[0].trim();
                String id = parts[1].trim();
                String smiles = parts[2].trim();
                int target = Integer.parseInt(parts[3].trim());
                
                double[] features = new double[17];
                for (int i = 0; i < 17; i++) {
                    features[i] = Double.parseDouble(parts[4 + i].trim());
                }
                
                samples.add(new Sample(ds, id, smiles, target, features));
                loadedLines++;
                
            } catch (NumberFormatException e) {
                System.err.println("\n=== ERRORE PARSING Riga " + lineNumber + " ===");
                System.err.println("Errore: " + e.getMessage());
                System.err.println("Riga: " + line);
                skippedLines++;
            }
        }
        br.close();
        
        System.out.println("Righe caricate: " + loadedLines);
        if (skippedLines > 0) {
            System.err.println("!!! ATTENZIONE: " + skippedLines + " righe saltate !!!");
        }
        System.out.println("======================\n");
        
        return samples;
    }
    
    public static void writeCSV(String filename, List<String> lines) throws IOException {
        PrintWriter pw = new PrintWriter(new FileWriter(filename));
        for (String line : lines) {
            pw.println(line);
        }
        pw.close();
    }
    
    public static String formatDouble(double value) {
        return String.format("%.6f", value).replace(",", ".");
    }
    
//    public static void main(String[] args) {
//        try {
//            System.out.println("=== GBT Molecular Classifier ===");
//            System.out.println("Parametri: trees=500, maxDepth=5, minNodeSize=3");
//            System.out.println("Features: 17 descriptors");
//
//            List<Sample> tsSamples = readCSV("TS.csv");
//            List<Sample> csSamples = readCSV("CS.csv");
//            List<Sample> vsSamples = readCSV("VS.csv");
//
//            System.out.printf("\nDataset caricati: TS=%d, CS=%d, VS=%d campioni%n",
//                            tsSamples.size(), csSamples.size(), vsSamples.size());
//
//            System.out.println("\nInizio training del modello...");
//            GradientBoostedTrees model = new GradientBoostedTrees();
//            model.train(tsSamples);
//
//            // EXPORT PMML
//            model.exportToPMML("model.pmml");
//
//            // TEST IMPORT PMML
//            System.out.println("\nTest di import PMML...");
//            GradientBoostedTrees loadedModel = GradientBoostedTrees.importFromPMML("model.pmml");
//
//            // Verifica che il modello caricato dia stessi risultati
//            System.out.println("Verifica predizioni modello caricato...");
//            double[] testFeatures = tsSamples.get(0).features;
//            double origProb = model.predictProbability(testFeatures);
//            double loadedProb = loadedModel.predictProbability(testFeatures);
//            System.out.printf("Probabilità originale: %.6f, Caricato: %.6f, Differenza: %.9f%n",
//                            origProb, loadedProb, Math.abs(origProb - loadedProb));
//
//            System.out.println("\nGenerazione predizioni...");
//
//            // Usa il modello originale per le predizioni
//            for (Sample s : tsSamples) {
//                s.probability = model.predictProbability(s.features);
//                s.prediction = model.predictClass(s.features);
//                s.confidence = Math.abs(s.probability - 0.5) * 2;
//            }
//
//            for (Sample s : csSamples) {
//                s.probability = model.predictProbability(s.features);
//                s.prediction = model.predictClass(s.features);
//                s.confidence = Math.abs(s.probability - 0.5) * 2;
//            }
//
//            for (Sample s : vsSamples) {
//                s.probability = model.predictProbability(s.features);
//                s.prediction = model.predictClass(s.features);
//                s.confidence = Math.abs(s.probability - 0.5) * 2;
//            }
//
//            List<String> predictionLines = new ArrayList<>();
//            predictionLines.add("DS,ID,SMILES,Target,Prediction,Probability,Confidence");
//
//            for (Sample s : tsSamples) {
//                predictionLines.add(String.format("%s,%s,%s,%d,%d,%s,%s",
//                                                s.ds, s.id, s.smiles, s.target,
//                                                s.prediction, formatDouble(s.probability), formatDouble(s.confidence)));
//            }
//
//            for (Sample s : csSamples) {
//                predictionLines.add(String.format("%s,%s,%s,%d,%d,%s,%s",
//                                                s.ds, s.id, s.smiles, s.target,
//                                                s.prediction, formatDouble(s.probability), formatDouble(s.confidence)));
//            }
//
//            for (Sample s : vsSamples) {
//                predictionLines.add(String.format("%s,%s,%s,%d,%d,%s,%s",
//                                                s.ds, s.id, s.smiles, s.target,
//                                                s.prediction, formatDouble(s.probability), formatDouble(s.confidence)));
//            }
//
//            writeCSV("predictions.csv", predictionLines);
//
//            ClassificationStats tsStats = new ClassificationStats(tsSamples);
//            ClassificationStats csStats = new ClassificationStats(csSamples);
//            ClassificationStats vsStats = new ClassificationStats(vsSamples);
//
//            System.out.println("\nRisultati:");
//            System.out.printf("TS - MCC: %s, Accuracy: %s%n", formatDouble(tsStats.mcc), formatDouble(tsStats.accuracy));
//            System.out.printf("CS - MCC: %s, Accuracy: %s%n", formatDouble(csStats.mcc), formatDouble(csStats.accuracy));
//            System.out.printf("VS - MCC: %s, Accuracy: %s%n", formatDouble(vsStats.mcc), formatDouble(vsStats.accuracy));
//
//            List<String> statsLines = new ArrayList<>();
//            statsLines.add("Dataset,BalancedAccuracy,Accuracy,Sensitivity,Specificity,PPV,NPV,MCC,TP,FP,TN,FN,OutOfDomain");
//
//            statsLines.add(String.format("TS,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,0",
//                                        formatDouble(tsStats.balancedAccuracy), formatDouble(tsStats.accuracy),
//                                        formatDouble(tsStats.sensitivity), formatDouble(tsStats.specificity),
//                                        formatDouble(tsStats.ppv), formatDouble(tsStats.npv), formatDouble(tsStats.mcc),
//                                        tsStats.tp, tsStats.fp, tsStats.tn, tsStats.fn));
//
//            statsLines.add(String.format("CS,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,0",
//                                        formatDouble(csStats.balancedAccuracy), formatDouble(csStats.accuracy),
//                                        formatDouble(csStats.sensitivity), formatDouble(csStats.specificity),
//                                        formatDouble(csStats.ppv), formatDouble(csStats.npv), formatDouble(csStats.mcc),
//                                        csStats.tp, csStats.fp, csStats.tn, csStats.fn));
//
//            statsLines.add(String.format("VS,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,0",
//                                        formatDouble(vsStats.balancedAccuracy), formatDouble(vsStats.accuracy),
//                                        formatDouble(vsStats.sensitivity), formatDouble(vsStats.specificity),
//                                        formatDouble(vsStats.ppv), formatDouble(vsStats.npv), formatDouble(vsStats.mcc),
//                                        vsStats.tp, vsStats.fp, vsStats.tn, vsStats.fn));
//
//            writeCSV("stats.csv", statsLines);
//
//            System.out.println("\n=== Completato! ===");
//            System.out.println("File generati: predictions.csv, stats.csv, model.pmml");
//
//        } catch (IOException e) {
//            System.err.println("Errore I/O: " + e.getMessage());
//            e.printStackTrace();
//        } catch (Exception e) {
//            System.err.println("Errore: " + e.getMessage());
//            e.printStackTrace();
//        }
//    }
}