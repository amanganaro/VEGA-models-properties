package insilico.devtox_pg.library;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author Alberto
 */
public class VirtualCompoundSubCategory implements Serializable{

    private int Index;
    private String SubIndex;
    private String Description;
    private ArrayList<VirtualCompound> Scaffolds;
    private ArrayList<VirtualCompound> Molecules;

    
    public VirtualCompoundSubCategory() {
        Index = 0;
        SubIndex = "";
        Description = "";
        Scaffolds = new ArrayList<>();
        Molecules = new ArrayList<>();        
    }
    

    /**
     * @return the Index
     */
    public int getIndex() {
        return Index;
    }

    /**
     * @param Index the Index to set
     */
    public void setIndex(int Index) {
        this.Index = Index;
    }

    /**
     * @return the SubIndex
     */
    public String getSubIndex() {
        return SubIndex;
    }

    /**
     * @param SubIndex the SubIndex to set
     */
    public void setSubIndex(String SubIndex) {
        this.SubIndex = SubIndex;
    }

    /**
     * @return the Description
     */
    public String getDescription() {
        return Description;
    }

    /**
     * @param Description the Description to set
     */
    public void setDescription(String Description) {
        this.Description = Description;
    }

    
    /**
     * @return the Scaffolds size
     */
    public int getScaffoldsSize() {
        return Scaffolds.size();
    }

    /**
     * @return the Scaffold
     */
    public ArrayList<VirtualCompound> getScaffolds() {
        return Scaffolds;
    }

    /**
     * @param index
     * @return the Scaffold no. index
     */
    public VirtualCompound getScaffold(int index) {
        return Scaffolds.get(index);
    }

    /**
     * @param VC
     */
    public void addScaffold(VirtualCompound VC) {
        Scaffolds.add(VC);
    }

    
    /**
     * @return the Molecules size
     */
    public int getMoleculesSize() {
        return Molecules.size();
    }

    /**
     * @return the Molecule
     */
    public ArrayList<VirtualCompound> getMolecules() {
        return Molecules;
    }

    /**
     * @param index
     * @return the Molecule no. index
     */
    public VirtualCompound getMolecule(int index) {
        return Molecules.get(index);
    }

    /**
     * @param VC
     */
    public void addMolecule(VirtualCompound VC) {
        Molecules.add(VC);
    }
    
}
