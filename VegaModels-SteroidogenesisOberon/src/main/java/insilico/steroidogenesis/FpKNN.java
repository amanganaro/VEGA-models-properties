package insilico.steroidogenesis;

import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.similarity.Tanimoto;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;

public class FpKNN {

    private final static int N_NEIGHBOURS = 2;
    private final static String KNN_DATA_URL = "/data/knn_training_smote.txt";

    private final Fingerprinter fprinter = new Fingerprinter();
    private KnnMol[] molecules;


    private class KnnMol {
        public String SMILES;
        public double Exp;
        public IBitFingerprint FP;
    }

    public FpKNN() throws InitFailureException {

        try {

            URL KnnDS = ismSteroidogenesis.class.getResource(KNN_DATA_URL);
            BufferedReader br = new BufferedReader(new InputStreamReader(KnnDS.openStream()));

            String line;
            int n = 0;
            while ((line = br.readLine()) != null)
                n++;
            br.close();

            molecules = new KnnMol[n];

            br = new BufferedReader(new InputStreamReader(KnnDS.openStream()));

            n = 0;
            while ((line = br.readLine()) != null) {
                String[] buf = line.split("\t");
                molecules[n] = new KnnMol();
                molecules[n].SMILES = buf[0];
                molecules[n].Exp = Double.valueOf(buf[1]);
                InsilicoMolecule mol = SmilesMolecule.Convert(molecules[n].SMILES);
                molecules[n].FP = fprinter.getBitFingerprint(mol.GetStructure());
                n++;
            }
            br.close();

        } catch (Exception e) {
            throw new InitFailureException(e);
        }

    }

    public double Calculate(InsilicoMolecule mol, boolean SkipExactMatch) throws Exception{

        IBitFingerprint fp = fprinter.getBitFingerprint(mol.GetStructure());

        double[] dist = new double[molecules.length];
        for (int i=0; i< molecules.length; i++) {
            dist[i] = 1.0 - Tanimoto.calculate(fp, molecules[i].FP);
        }

        int[] neigh = getSortedNeighbours(dist,N_NEIGHBOURS, SkipExactMatch);

        double pred = -1;
        if (molecules[neigh[0]].Exp == molecules[neigh[1]].Exp) {
            pred = molecules[neigh[0]].Exp;
        } else {
            // discordance - at least one is positive, and we use this as prediction
            // but if delta distances is > 0.05 we use just the first one
            double delta = Math.abs(dist[neigh[0]] - dist[neigh[1]]);
            if (delta > 0.05)
                pred = molecules[neigh[0]].Exp;
            else
                pred = 1;
        }

        return pred;
    }


    private int[] getSortedNeighbours(double[] DistanceArray, int nNeigh, boolean SkipExactMatch) {

        double[] distances = DistanceArray.clone();
        int[] Res = new int[nNeigh];

        for (int neigh = 0; neigh < nNeigh; neigh++) {

            int minIndex = 0;
            if (SkipExactMatch)
                if (distances[0] == 0)
                    minIndex = 1;

            for (int i=0; i< distances.length; i++) {
                if (SkipExactMatch)
                    if (distances[i] == 0)
                        continue;
                if (distances[i] < distances[minIndex])
                    minIndex = i;
            }
            Res[neigh] = minIndex;
            distances[minIndex] = Double.POSITIVE_INFINITY;
        }

        return Res;
    }
}
