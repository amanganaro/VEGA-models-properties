package insilico.devtox_pg.library;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author Alberto
 */
public class VirtualCompoundCategory implements Serializable {

    private int Index;
    private String Description;
    private ArrayList<VirtualCompoundSubCategory> SubCategory;

    public VirtualCompoundCategory() {
        Index = 0;
        Description = "";
        SubCategory = new ArrayList<>();
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
     * @return the SubCategory
     */
    public ArrayList<VirtualCompoundSubCategory> getSubCategory() {
        return SubCategory;
    }

    public void addSubCategory(VirtualCompoundSubCategory cat) {
        this.SubCategory.add(cat);
    }     
}
