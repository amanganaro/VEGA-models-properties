package insilico.bcf_caesar.descriptors;

import java.io.*;
import java.util.ArrayList;
import java.util.List;


/**
 * Utility class used to read csv dataset from file
 * 
 * @version 1.0
 * @author  Alberto Manganaro,
 *          Laboratory of Environmental Toxicology and Chemistry,
 *          Istituto di Ricerche Farmacologiche "Mario Negri", Milano, Italy
 */
public class RBFNNDataReader {

    private class StringSplitter {

        private int nfield; 
        private char fieldSep;
        private List<String> field = new ArrayList<String>();
        private String fld = new String();


        /**
         * Constructor, field separator is set to TAB
         */
        public StringSplitter() {
            fieldSep = '\t';
        }


        /**
         * Constructor 
         * 
         * @param sep   character to be use as field separator
         */
        public StringSplitter(char sep) {
            fieldSep = sep;
        }


        /**
         * Parses the given line and returns the number of tokens
         * 
         * @param line  string to be parsed
         * @return      the number of tokens read 
         */
        public int split(String line) {

            int i, j = 0;
            fld = "";

            field.clear();

            nfield = 0;

            if (line.length() == 0) 
                return 0;

            i = 0;

            do {
                if (i < line.length() )
                    j = advplain(line, i); 

                field.add(fld); 
                nfield++; 
                i = j + 1; 
            } while (j < line.length()); 

            return nfield; 

        }


        /**
         * Gives the n-th token of the parsed string
         * 
         * @param n     number of the token to be returned
         * @return      n-th token
         */
        public String getfield(int n) { //This function will return the field values

            if (n < 0 || n >= nfield) 
                return ""; 
            else 
                return field.get(n).toString(); 

        }


        private int advplain(String line, int i) {

            int j;

            j = line.indexOf(fieldSep, i);

            if (j == -1) {
                fld = line.substring(i); 
                return line.length(); 
            } else
                fld = line.substring(i,j); 

            return j; 

        } 

    }
    
    
    public double[][] Data;
    public String[] Id;
    public int dim_i;
    public int dim_j;


    /**
     * Constructor to be used with a stream input
     * 
     * @param Filename
     * @param SkipLine
     * @param SkipColumns
     * @param IdColumn
     * @throws IOException
     */
    public RBFNNDataReader(InputStream Filename, int SkipLine, int SkipColumns,
                           int IdColumn) throws IOException {

        DataInputStream in = new DataInputStream(Filename);
        ExecReader(in, SkipLine, SkipColumns, IdColumn);
        
    }

    
    /**
     * Constructor to be used with a file input
     * 
     * @param Filename
     * @param SkipLine
     * @param SkipColumns
     * @param IdColumn
     * @throws IOException
     */
    public RBFNNDataReader(String Filename, int SkipLine, int SkipColumns,
                           int IdColumn) throws IOException {

        FileInputStream fstream = new FileInputStream(Filename);
        DataInputStream in = new DataInputStream(fstream);
        ExecReader(in, SkipLine, SkipColumns, IdColumn);
        
    }
    
    
    private void ExecReader(DataInputStream in, int SkipLine, int SkipColumns,
                            int IdColumn) throws IOException {

	StringSplitter parser = new StringSplitter(',');
	int n=0, i=0, j=0, k=0;
        double[][] BufArray = new double[1][1];
        String[] BufVector = new String[1];
        Id = new String[1];

        dim_i = 0;
        if (SkipColumns<0)
            SkipColumns = 0;

        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;

            // Skips required lines
            if (SkipLine>0)
            	for (i=0; i<SkipLine; i++)
            		strLine = br.readLine();
            
            while ((strLine = br.readLine()) != null)   {

                // Gets number of tokens
                n = parser.split(strLine);

                // Resizes buffer array for the new line to be added
               	BufArray = new double[dim_i+1][n-SkipColumns];
               	for (i=0;i<dim_i;i++)
               		for (j=0;j<(n-SkipColumns);j++)
               			BufArray[i][j] = Data[i][j];

                // Adds the new line read
               	for (j=0;j<n-SkipColumns;j++) {
               		BufArray[dim_i][j] = Double.parseDouble(parser.getfield(j + SkipColumns));
               	}

                // Resizes buffer vector for the new id field
               	BufVector = new String[dim_i+1];
               	for (i=0;i<dim_i;i++)
           			BufVector[i] = Id[i];

                // Adds the new id field
                if (IdColumn>0)
               		BufVector[dim_i] = parser.getfield(IdColumn-1);
                else
               		BufVector[dim_i] = String.valueOf(dim_i+1).toString();

                Data = BufArray;
                Id = BufVector;
                dim_i++;

            }

            dim_j = n-SkipColumns;
            in.close();

        } catch (Exception e){
            throw new IOException("Unable to read from file");
        }

    }
    
    
    public int[][] DataAsInt() {
    	
    	int[][] BufData = new int[dim_i][dim_j];
    	
        for (int i=0; i<dim_i; i++)
        	for (int j=0; j<dim_j; j++)
        		BufData[i][j] = (int) Data[i][j];
    	
        return BufData;
    }

}
