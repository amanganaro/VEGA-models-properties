package insilico.steroidogenesis;

import insilico.core.exception.InitFailureException;
import org.junit.runners.model.InitializationError;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Arrays;

public class KNN {

    private final double[][] dataX;
    private final int[] dataY;

    double[] sample;


    public KNN(URL KnnData) throws InitFailureException {

        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(KnnData.openStream()));

            String line = br.readLine(); // headers
            int col = line.split("\t").length - 1;

            int n = 0;
            while ((line = br.readLine()) != null)
                n++;
            br.close();

            this.dataX = new double[n][col];
            this.dataY = new int[n];

            br = new BufferedReader(new InputStreamReader(KnnData.openStream()));

            line = br.readLine(); // headers
            n = 0;
            while ((line = br.readLine()) != null) {
                String[] buf = line.split("\t");
                this.dataY[n] = Integer.valueOf(buf[0]);
                for (int c = 0; c < col; c++)
                    this.dataX[n][c] = Double.valueOf(buf[c + 1]);
                n++;
            }
        } catch (Throwable e) {
            throw new InitFailureException("Unable to init KNN data - " + e.getMessage());
        }

    }


    public int getPrediction(double[] sample) {
        return getPrediction(sample, null);
    }


    public int getPrediction(double[] sample, Integer SkipSample){
        this.sample = sample;
        double[] distances = calculateDistances(sample);
        int[] bestNeigh = getSortedNeighbours(distances, 2, SkipSample);

        // matching exp values
        if (dataY[bestNeigh[0]] == dataY[bestNeigh[1]]) {
            System.out.println(dataY[bestNeigh[0]] + "\t" + "matching");
            return dataY[bestNeigh[0]];
        }

        // otherwise, not matching but both with same distance
        if (distances[bestNeigh[0]] == distances[bestNeigh[1]]) {
            // always wins inactive prediction
            System.out.println(dataY[bestNeigh[1]] + "\t" + "eq dist");
            return dataY[bestNeigh[1]];
//            return 0;
        }

        // otherwise, prediction come from first (less distant) neigh
        System.out.println(dataY[bestNeigh[0]] + "\t" + "best neigh");
        return dataY[bestNeigh[0]];

//        // if first neigh has distance = 0, take its value
//        if (distances[bestNeigh[0]] == 0)
//            return dataY[bestNeigh[0]];
//
//        // otherwise weights experimentals by distance
//        double[] ClassWeights = new double[2];
//        ClassWeights[0] = 1000; ClassWeights[1] = 1000;
//        for (int simIndex : bestNeigh)
//            ClassWeights[dataY[simIndex]] -= distances[simIndex];
//
//        if (ClassWeights[0] < ClassWeights[1])
//            return 0;
//        else
//            return 1;

    }


    public int getPredictionCustom(double[] sample, Integer SkipSample){
        this.sample = sample;
        double[] distances = calculateDistances(sample);
        int[] bestNeigh = getSortedNeighbours(distances, 3, SkipSample);

        // matching exp values
        if ( (dataY[bestNeigh[0]] == dataY[bestNeigh[1]]) && (dataY[bestNeigh[0]] == dataY[bestNeigh[2]])) {
            System.out.println(dataY[bestNeigh[0]] + "\t" + "matching");
            return dataY[bestNeigh[0]];
        }


        // otherwise weights experimentals by distance
        double[] ClassWeights = new double[2];
        ClassWeights[0] = 0; ClassWeights[1] = 0;
        for (int simIndex : bestNeigh)
            ClassWeights[dataY[simIndex]] += distances[simIndex];

        if (ClassWeights[0] < ClassWeights[1]) {
            System.out.println("0");
            return 0;
        } else {
            System.out.println("1");
            return 1;
        }
    }

    private int[] getSortedNeighbours(double[] DistanceArray, int nNeigh, Integer SkipSample) {

        double[] distances = DistanceArray.clone();
        int[] Res = new int[nNeigh];

        for (int neigh = 0; neigh < nNeigh; neigh++) {
            int minIndex = 0;
            if (SkipSample != null)
                if (SkipSample == 0)
                    minIndex = 1;
            for (int i=0; i< distances.length; i++) {
                if (SkipSample != null)
                    if (SkipSample == i)
                        continue;
                if (distances[i] < distances[minIndex])
                    minIndex = i;
            }
            Res[neigh] = minIndex;
            distances[minIndex] = Double.POSITIVE_INFINITY;
        }

        return Res;
    }

    private double[] calculateDistances(double[] sample) {
        double[] distances = new double[dataX.length];
        for (int i = 0; i < dataX.length; i++) {
            distances[i] = calculateDistance(dataX[i], sample);
        }
        return distances;
    }

    private double calculateDistance(double[] point, double[] sample) {
        double dist = 0;
        for (int i = 0; i < point.length; i++)
            dist += Math.pow(point[i] - sample[i], 2);
        return Math.sqrt(dist);
    }

}
