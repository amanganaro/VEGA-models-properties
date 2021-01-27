package insilico.devtox_pg.library;

import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.MDLMolecule;
import insilico.core.molecule.conversion.file.MoleculeFile;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Reader for SDF file. The reader can recognize tags for molecule's id and
 * CAS number, if provided.
 * 
 * MODIFIED for PG dataset where some molecules have bond order = 5
 * which means a single or a double bond
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class PGMoleculeFileSDF extends MoleculeFile {

    private String CASTag;
    private String IdTag;

    
    /**
     * Constructor.
     */
    public PGMoleculeFileSDF() {
        super();
        CASTag = null;
        IdTag = null;
    }
    
    
    /**
     * Reads and return next molecule found in file. Returns an IOException if
     * the file has not been already opened. It always return the result of 
     * parsed lines, if they contain invalid SMILES this will be eventually
     * seen with the isValid() method of the returned molecule.
     * 
     * @return the next parsed molecule in the file
     * @throws IOException
     */
    @Override
    public InsilicoMolecule ReadNext() throws IOException {
        
        if (!isFileOpen)
            throw new IOException("File is not open");
        
        try {
            // Parses SDF and separates each single MDL molecule
            String MDLMol = "";
            boolean Proceed = true;

            while (Proceed) {
                String CurLine = reader.readLine();

                if ((CurLine == null) || (CurLine.compareTo("$$$$") == 0)) {
                    Proceed = false;
                } else {
                    MDLMol += CurLine + "\n";
                }
            }

            if (MDLMol.isEmpty())
                return null;

            Count++;

            byte[] bytes = MDLMol.getBytes();
            InsilicoMolecule m = MDLMolecule.Convert(bytes, CASTag, IdTag);
            if (m.GetId().isEmpty())
                m.SetId("Molecule " + Count);
            return m;
        } catch (IOException e) {
            log.error("Error while reading file " + this.FileName + " (" + e.getMessage() + ")");
            throw(e);
        }

    }
    
    
    /**
    * 
    */
    public ArrayList<InsilicoMolecule> ReadNextMultiple() throws IOException {
        
        if (!isFileOpen)
            throw new IOException("File is not open");
        
        try {
            // Parses SDF and separates each single MDL molecule
            String MDLMol = "";
            boolean Proceed = true;

            while (Proceed) {
                String CurLine = reader.readLine();

                if ((CurLine == null) || (CurLine.compareTo("$$$$") == 0)) {
                    Proceed = false;
                } else {
                    MDLMol += CurLine + "\n";
                }
            }

            if (MDLMol.isEmpty())
                return null;

            ArrayList<InsilicoMolecule> res = new ArrayList<>();
            
            if (CheckForErrBond(MDLMol)) {
                String MDL_1 = SetBonds(MDLMol, 5, 1);
                String MDL_2 = SetBonds(MDLMol, 5, 2);
                
                byte[] bytes_1 = MDL_1.getBytes();
                InsilicoMolecule m_1 = MDLMolecule.Convert(bytes_1, CASTag, IdTag);
                if (m_1.GetId().isEmpty())
                    m_1.SetId("Molecule " + Count + "_bnd1");
                
                byte[] bytes_2 = MDL_2.getBytes();
                InsilicoMolecule m_2 = MDLMolecule.Convert(bytes_2, CASTag, IdTag);
                if (m_2.GetId().isEmpty())
                    m_2.SetId("Molecule " + Count + "_bnd2");
                
                res.add(m_1);
                res.add(m_2);
            } else {
                byte[] bytes = MDLMol.getBytes();
                InsilicoMolecule m = MDLMolecule.Convert(bytes, CASTag, IdTag);
                if (m.GetId().isEmpty())
                    m.SetId("Molecule " + Count);
                res.add(m);    
            }
            
            Count++;

            return res;
            
        } catch (IOException e) {
            log.error("Error while reading file " + this.FileName + " (" + e.getMessage() + ")");
            throw(e);
        }

    }

    
    private boolean CheckForErrBond(String MDL) {
        String[] MDLlines = MDL.split("\n");

        int idx = 4; // start after header lines

        // try to skip lines until the bond section
        while(idx < MDLlines.length) {
            String[] line = MDLlines[idx].trim().split(" +"); // split on multiple spaces
            if (line.length<4) break;
            if (isNumeric(line[3])) break;
            idx++;
        }

        // check bonds
        while(idx < MDLlines.length) {
            String[] line = MDLlines[idx].trim().split(" +"); // split on multiple spaces
            if (line.length<3) break;
            int bond = Integer.parseInt(line[2]);
            if (bond == 5)
                return true;
            idx++;
        }

        return false;
    }
    
    private String SetBonds(String MDL, int sourceBond, int newBond) {
        
        String[] MDLlines = MDL.split("\n");
        ArrayList<String> newMDL = new ArrayList<>();

        // headers
        int idx = 0;
        while(idx<4) {
            newMDL.add(MDLlines[idx]);
            idx++;
        }
        
        // try to skip lines until the bond section
        while(idx < MDLlines.length) {
            String[] line = MDLlines[idx].trim().split(" +"); // split on multiple spaces
            if (line.length<4) break;
            if (isNumeric(line[3])) break;
            newMDL.add(MDLlines[idx]);
            idx++;
        }

        // check bonds
        while(idx < MDLlines.length) {
            String curLine = MDLlines[idx];
            String[] line = MDLlines[idx].trim().split(" +"); // split on multiple spaces
            if (line.length>=3) {
                int bond = Integer.parseInt(line[2]);
                if (bond == sourceBond) {
                    int chrIdx = 0;
                    boolean intoChar = false;
                    int colIdx = 0;
                    while(chrIdx < curLine.length()) {
                        if (curLine.charAt(chrIdx) == ' ') {
                            if (intoChar) {
                                intoChar = false;
                            }
                        } else {
                            if (!intoChar) {
                                intoChar = true;
                                colIdx++;
                            }
                        }
                        
                        if (colIdx == 3) break;
                        
                        chrIdx++;
                    }
                    
                    char[] buf = curLine.toCharArray();
                    buf[chrIdx] = String.valueOf(newBond).charAt(0);
                    curLine = String.valueOf(buf);                    

                }
            }
            newMDL.add(curLine);
            idx++;
        }

        String res = "";
        for (String s : newMDL)
            res += s + "\n";
        return res;
    }    
    
    public static boolean isNumeric(String str) {  
        try {  
          double d = Double.parseDouble(str);  
        } catch(NumberFormatException nfe) {  
          return false;  
        }  
        return true;  
    }    
    
    
    /**
     * Reads all the molecules in the file and returns them as an ArrayList. 
     * Note: this method will start reading from the current file mark, thus
     * if some molecules have been already parsed with the ReadNext() method, 
     * this method will return all the molecules starting from that point.
     * 
     * @return array of parsed molecules
     * @throws IOException 
     */
    @Override
    public ArrayList<InsilicoMolecule> ReadAll() throws IOException {

        ArrayList<InsilicoMolecule> Mols = new ArrayList<>();
        ArrayList<InsilicoMolecule> curmols;
        while ((curmols=this.ReadNextMultiple()) != null) {
            for (InsilicoMolecule m : curmols)
                Mols.add(m);
        }
        return Mols;
            
    }

    /**
     * @param CASTag the CASTag to set
     */
    public void setCASTag(String CASTag) {
        this.CASTag = CASTag;
    }

    /**
     * @param IdTag the IdTag to set
     */
    public void setIdTag(String IdTag) {
        this.IdTag = IdTag;
    }

    
    
}
