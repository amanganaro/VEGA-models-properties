package insilico.fish_nic_tk.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;

/**
 * Topological distances descriptors.
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class TopologicalDistances extends DescriptorBlock {

    private final static long serialVersionUID = 1L;

    Logger logger = LoggerFactory.getLogger(TopologicalDistances.class);

    private final static String BlockName = "Topological Distances Descriptors";
    
    private final static int MAX_TOPO_DISTANCE = 10;

    private final static String AtomCouples[][] = { 
        {"C", "C"},
        {"C", "N"},
        {"C", "O"},
        {"C", "S"},
        {"C", "P"},
        {"C", "F"},
        {"C", "Cl"},
        {"C", "Br"},
        {"C", "I"},
        {"N", "N"},
        {"N", "O"},
        {"N", "S"},
        {"N", "P"},
        {"N", "F"},
        {"N", "Cl"},
        {"N", "Br"},
        {"N", "I"},
        {"O", "O"},
        {"O", "S"},
        {"O", "P"},
        {"O", "F"},
        {"O", "Cl"},
        {"O", "Br"},
        {"O", "I"},
        {"S", "S"},
        {"S", "P"},
        {"S", "F"},
        {"S", "Cl"},
        {"S", "Br"},
        {"S", "I"},
        {"P", "P"},
        {"P", "F"},
        {"P", "Cl"},
        {"P", "Br"},
        {"P", "I"},
        {"F", "F"},
        {"F", "Cl"},
        {"F", "Br"},
        {"F", "I"},
        {"Cl", "Cl"},
        {"Cl", "Br"},
        {"Cl", "I"},
        {"Br", "Br"},
        {"Br", "I"},
        {"I", "I"}
    };
    
    
    
    /**
     * Constructor. This should not be used, no weight is specified. The 
     * overloaded constructors should be used instead.
     */    
    public TopologicalDistances() {
        super();
        this.Name = TopologicalDistances.BlockName;
    }

    
    
    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        for (int i=0; i<AtomCouples.length; i++)
            Add("T(" + AtomCouples[i][0] + ".." + AtomCouples[i][1] + ")", "");
        for (int lag=1; lag<= MAX_TOPO_DISTANCE; lag++)
            for (int i=0; i<AtomCouples.length; i++) {
                Add("B"+lag+"(" + AtomCouples[i][0] + ".." + AtomCouples[i][1] + ")", "Presence/absence of "+ AtomCouples[i][0] + "-" + AtomCouples[i][1] + " at topological distance " + lag);
                Add("F"+lag+"(" + AtomCouples[i][0] + ".." + AtomCouples[i][1] + ")", "Frequency of "+ AtomCouples[i][0] + "-" + AtomCouples[i][1] + " at topological distance " + lag);
            }
        SetAllValues(Descriptor.MISSING_VALUE);
    }
    

        
    /**
     * Calculate descriptors for the given molecule.
     * 
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {

        // Generate/clears descriptors
        GenerateDescriptors();

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            logger.warn(e.getMessage());
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }
        
        for (int d=0; d<AtomCouples.length; d++) {
        
            int descT = 0;
            int[] descB = new int[MAX_TOPO_DISTANCE];
            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descB, 0);
            Arrays.fill(descF, 0);
            
            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouples[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouples[d][1])) {

                            // T (sum of topo distances)
                            if (TopoMat[i][j] > 2) // DA VEDERE PERCHE MAGGIORE DI 2
                                    descT += TopoMat[i][j];
                            
                            // B (presence of pair) and F (number of couples)
                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descB[TopoMat[i][j]-1] = 1;
                                descF[TopoMat[i][j]-1]++;
                            }
                            
                        }                        
                    }
                }
            }      
            
            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouples[d][0].compareTo(AtomCouples[d][1]) == 0) {
                descT /= 2;
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }
            
            SetByName("T(" + AtomCouples[d][0] + ".." + AtomCouples[d][1] + ")", descT);
            for (int i=0; i<descB.length; i++)            
                SetByName("B"+(i+1)+"(" + AtomCouples[d][0] + ".." + AtomCouples[d][1] + ")", descB[i]);
            for (int i=0; i<descF.length; i++)            
                SetByName("F"+(i+1)+"(" + AtomCouples[d][0] + ".." + AtomCouples[d][1] + ")", descF[i]);
                
        }

    }


    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException 
     */
    @Override
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        TopologicalDistances block = new TopologicalDistances();
        block.CloneDetailsFrom(this);
        return block;
    }

    
}
