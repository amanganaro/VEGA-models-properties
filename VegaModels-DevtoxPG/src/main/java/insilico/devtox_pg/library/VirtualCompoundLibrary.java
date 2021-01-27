package insilico.devtox_pg.library;

import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
//import insilico.core.tools.logger.InsilicoLogger;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;

/**
 *
 * @author Alberto
 */
public class VirtualCompoundLibrary {

    //
    // Process a given directory with SDF files organized as provided
    // by PG (multiple directories with categories)
    // and saves all molecules in SMI files (in the given dest directory)
    //
    public static void ProcessDirectorySDF(Path SourceDir, String DestDir) throws InitFailureException {

        try(DirectoryStream<Path> stream = Files.newDirectoryStream(SourceDir)) {
            for (Path path : stream) {
                if(path.toFile().isDirectory()) {
                    ProcessDirectorySDF(path, DestDir);
                } else {

                    System.out.println("* Processing: " + path.toAbsolutePath().toString());
                    String sdfname = path.getFileName().toString();
                    String sminame = (sdfname.split("\\."))[0] + ".smi";
                    
                    PrintWriter f = new PrintWriter(new File(DestDir + "\\" + sminame));
                    
                    PGMoleculeFileSDF sdf = new PGMoleculeFileSDF();
                    sdf.OpenFile(path.toAbsolutePath().toString());
                    
                    ArrayList<InsilicoMolecule> mols = sdf.ReadAll();

                    for (int i=0; i<mols.size(); i++) {
                        InsilicoMolecule m = mols.get(i);
                        if (!m.IsValid()) 
                            System.out.println((i+1) + " INVALID molecule - " + m.GetErrors().GetMessages());
                        else
                            f.println(m.GetSMILES());
                    }
                    
                    f.close();
                }
            }
        } catch(IOException e) {
            e.printStackTrace();
        }        
    }
    
    
    //
    // Process a directory with all SMI files as created by the previous method
    // and creates the VirtualCompoundCategory objects with parsed molecules
    //
    public static void ProcessDirectorySMI(Path SourceDir) throws InitFailureException, FileNotFoundException, IOException {

        HashMap<String,VirtualCompoundSubCategory> SubCategories = new HashMap<>();
        
        try(DirectoryStream<Path> stream = Files.newDirectoryStream(SourceDir)) {
            for (Path path : stream) {
                if(!path.toFile().isDirectory()) {
                    String name = path.getFileName().toString();
                    String cat;
                    if (name.split(" ").length > 1) 
                        cat = name.split(" ")[0];
                    else
                        cat = name.split("\\.")[0];
                    
                    // get number and sub-index
                    int idx = 0;
                    while (idx < cat.length()) {
                        if ( ! Character.isDigit(cat.charAt(idx)) )
                            break;
                        idx++;
                    }
                    int curIndex = Integer.parseInt(cat.substring(0, idx));
                    String curSubIndex = cat.substring(idx);
                    
                    // check if new category should be created
                    if (!SubCategories.containsKey(cat)) {
                        VirtualCompoundSubCategory VCC = new VirtualCompoundSubCategory();
                        VCC.setIndex(curIndex);
                        VCC.setSubIndex(curSubIndex);
                        SubCategories.put(cat, VCC);
                    }

                    // read current file
                    BufferedReader br = new BufferedReader(new FileReader(path.toFile()));
                    String line;
                    ArrayList<String> Mols = new ArrayList<>();
                    while ( (line = br.readLine())!= null )
                        if (line.length()>0)
                        Mols.add(line);
                    
                    // add mols or scaffold to category
                    VirtualCompoundSubCategory VCC = SubCategories.get(cat);
                    for (String s : Mols) {
                        QueryAtomContainer Query;
                        try {
                            if ( (name.contains("Scaffold")) || (name.contains("scaffold"))) {
                                VirtualCompound VC = new VirtualCompound(s, true);
                                VCC.addScaffold(VC);
                            } else {
                                VirtualCompound VC = new VirtualCompound(s, false);
                                VCC.addMolecule(VC);
                            }
                        } catch (Throwable ex) {
                            System.out.println("Invalid molecule (as SMARTS) in " + cat + ": " + s + "(" + ex.getMessage() + ")");
                        }
                    }
                    SubCategories.put(cat, VCC);
                }
            }
        } catch(IOException e) {
            e.printStackTrace();
        }  
        
        // organize into categories
        HashMap<Integer,VirtualCompoundCategory> Categories = new HashMap<>();        
        for (String s : SubCategories.keySet()) {
            VirtualCompoundSubCategory VCC = SubCategories.get(s);
            int catIndex = VCC.getIndex();
            VirtualCompoundCategory curCat;
            if (Categories.containsKey(catIndex))
                curCat = Categories.get(catIndex);
            else {
                curCat = new VirtualCompoundCategory();
                curCat.setIndex(catIndex);
            }
            curCat.addSubCategory(VCC);
            Categories.put(catIndex, curCat);
        }
        
        // serialize
        for (int s : Categories.keySet()) {
            VirtualCompoundCategory VCC = Categories.get(s);
            FileOutputStream fos = new FileOutputStream(VCC.getIndex() + ".dat" );
            ObjectOutputStream out = new ObjectOutputStream(fos);
            out.writeObject(VCC);
            out.close();
        }        
        
        // Print final stats
        System.out.println();
        System.out.println();
        for (Integer i : Categories.keySet()) {
            VirtualCompoundCategory VCC = Categories.get(i);
            System.out.println("----------------------------------------");
            System.out.println("Category " + i + ": ");
            System.out.println("Index\tSubindex\tScaffolds\tMolecules");
            for (VirtualCompoundSubCategory SVCC : VCC.getSubCategory()) {
                System.out.println(VCC.getIndex() + "\t" + SVCC.getSubIndex() + "\t" + SVCC.getScaffoldsSize() + "\t" + SVCC.getMoleculesSize() + "\t");
            }
        }
    }
    
    
}
