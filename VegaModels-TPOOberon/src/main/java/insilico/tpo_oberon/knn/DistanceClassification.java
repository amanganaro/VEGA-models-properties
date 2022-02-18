package insilico.tpo_oberon.knn;

public class DistanceClassification {

    private double distance;
    private int classification;
    private double weight;

    public DistanceClassification(double distance, int classification) {
        this.distance = distance;
        this.classification = classification;
        weight = 1 / distance;
    }

    public double getDistance() {
        return distance;
    }

    public void setDistance(double distance) {
        this.distance = distance;
    }

    public int getClassification() {
        return classification;
    }

    public void setClassification(int classification) {
        this.classification = classification;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }
}
