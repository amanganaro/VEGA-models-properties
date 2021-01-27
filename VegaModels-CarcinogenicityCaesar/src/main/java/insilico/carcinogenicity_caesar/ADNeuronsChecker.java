package insilico.carcinogenicity_caesar;

import insilico.core.exception.InitFailureException;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;

/**
 *
 * @author amanganaro
 */
public class ADNeuronsChecker {

    public int[][] NeuronsPos;
    public int[][] NeuronsNeg;

    public int Population;
    public double Concordance;


    public ADNeuronsChecker(int nX, int nY, String FileName) throws InitFailureException {

        NeuronsPos = new int[nX][nY];
        NeuronsNeg = new int[nX][nY];
        for (int i=0; i<nX; i++) 
            for (int j=0; j<nY; j++) {
                NeuronsNeg[i][j] = 0;
                NeuronsPos[i][j] = 0;
            }

        try {

            URL uNet = getClass().getResource(FileName);
            DataInputStream in = new DataInputStream(uNet.openStream());
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
//            StringSplitter parser = new StringSplitter('\t');
            String strLine;

            br.readLine();

            while ((strLine = br.readLine())!=null) {
                String[] buf = strLine.split("\t");
//                parser.split(strLine);
                if ( (buf[0].compareTo("-")==0) ||
                     (buf[1].compareTo("-")==0) ||
                     (buf[2].compareTo("-")==0) )
                    continue; // compound not calculated
                int CurPr = Integer.parseInt(buf[0]);
                int X = Integer.parseInt(buf[1]);
                int Y = Integer.parseInt(buf[2]);
                if (CurPr == 1) 
                    NeuronsPos[X][Y]++;
                else
                    NeuronsNeg[X][Y]++;
            }
            in.close();          

        } catch (IOException | NumberFormatException e) {
            throw new InitFailureException("Unable to init neurons data for AD check.");
        }    

        Population = 0;
        Concordance = 0;
    }


    public void Check(int nX, int nY, int PredictedClass) {

        Population = NeuronsPos[nX][nY] + NeuronsNeg[nX][nY];
        Concordance = 0;
        if (Population>0) {
            if (PredictedClass == 1) 
                Concordance = ((double)NeuronsPos[nX][nY]) / ((double)Population);
            else
                Concordance = ((double)NeuronsNeg[nX][nY]) / ((double)Population);
        }            
    }
    
}
