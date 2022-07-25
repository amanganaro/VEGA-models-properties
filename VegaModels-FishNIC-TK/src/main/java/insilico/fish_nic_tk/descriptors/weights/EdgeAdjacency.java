package insilico.fish_nic_tk.descriptors.weights;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Edge Adjacency descriptors.
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class EdgeAdjacency extends DescriptorBlock {

    private final static long serialVersionUID = 1L;
    private final static String BlockName = "Edge Adjacency Descriptors";

    private final static int MinEig = 1;
    private final static int MaxEig = 15;
    private final static int MinSM = 1;
    private final static int MaxSM = 15;
    
    public final static String PARAMETER_WEIGHT_X = "weightX";  // Edge degree
    public final static String PARAMETER_WEIGHT_D = "weightD";  // Dipole moments
    public final static String PARAMETER_WEIGHT_R = "weightR";  // Resonance integral
    
    private final static short WEIGHT_X_IDX = 0;
    private final static short WEIGHT_D_IDX = 1;
    private final static short WEIGHT_R_IDX = 2;
    private final static String[] WEIGHT_SYMBOL = {"ed", "dm", "ri"};

    private boolean defaultDescriptors;



    /**
     * Constructor. This should not be used, no weight is specified. The 
     * overloaded constructors should be used instead.
     */
    public EdgeAdjacency() {
        super();
        this.Name = EdgeAdjacency.BlockName;
        this.defaultDescriptors = true;
    }

    public EdgeAdjacency(boolean defaultDescriptors) {
        super();
        this.Name = EdgeAdjacency.BlockName;
        this.defaultDescriptors = defaultDescriptors;
    }

    
    
    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        ArrayList<Integer> weightList = BuildWeightList();
        for (Integer curWeight : weightList) {
            Add("SpMax" + WEIGHT_SYMBOL[curWeight], "");
            for (int i=MinEig; i<=MaxEig; i++) 
                Add("EEig" + i + WEIGHT_SYMBOL[curWeight], "");
            for (int i=MinSM; i<=MaxSM; i++) {
                Add("ESpm" + i + WEIGHT_SYMBOL[curWeight], "");
                Add("ESpm" + i + WEIGHT_SYMBOL[curWeight] + "D5", ""); // compatibility with Dragon 5
            }
        }
        SetAllValues(Descriptor.MISSING_VALUE);
    }


    private ArrayList<Integer> BuildWeightList() {
        ArrayList<Integer> w = new ArrayList<>();
        if (defaultDescriptors) {
            w.add((int) WEIGHT_X_IDX);
            w.add((int) WEIGHT_D_IDX);
            w.add((int) WEIGHT_R_IDX);
        } else {
            if (getBoolProperty(PARAMETER_WEIGHT_X))
                w.add((int) WEIGHT_X_IDX);
            if (getBoolProperty(PARAMETER_WEIGHT_D))
                w.add((int) WEIGHT_D_IDX);
            if (getBoolProperty(PARAMETER_WEIGHT_R))
                w.add((int) WEIGHT_R_IDX);
        }

        return w;
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
        
        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }
        
        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }
        
        Matrix DataMatrix = null;

        // Generates weight
        ArrayList<Integer> weightList = BuildWeightList();

        // Cycle for all weights
        for (Integer curWeight : weightList) {
           
            if (curWeight == WEIGHT_X_IDX) {
                double[][] EdgeDegreeMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
                for (int i=0; i<EdgeAdjMat.length; i++)
                    for (int j=0; j<EdgeAdjMat[0].length; j++)
                        EdgeDegreeMat[i][j] = EdgeAdjMat[i][j][1];

                DataMatrix = new Matrix(EdgeDegreeMat);
            } 
            
            if (curWeight == WEIGHT_D_IDX) {
                double[][] EdgeDipoleMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
                for (int i=0; i<EdgeAdjMat.length; i++)
                    for (int j=0; j<EdgeAdjMat[0].length; j++)
                        EdgeDipoleMat[i][j] = EdgeAdjMat[i][j][0];

                for (int i=0; i<m.getBondCount(); i++) {
                    IAtom a =  m.getBond(i).getAtom(0);
                    IAtom b =  m.getBond(i).getAtom(1);

                    double CurVal = GetDipoleMoment(m, a, b);
                    if (CurVal == 0)
                        CurVal = GetDipoleMoment(m, b, a);
                    EdgeDipoleMat[i][i] = CurVal;
                }

                DataMatrix = new Matrix(EdgeDipoleMat);
            } 
            
            if (curWeight == WEIGHT_R_IDX) {
                double[][] EdgeResMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
                for (int i=0; i<EdgeResMat.length; i++)
                    for (int j=0; j<EdgeResMat[0].length; j++)
                        EdgeResMat[i][j] = EdgeAdjMat[i][j][1];

                for (int i=0; i<m.getBondCount(); i++) {
                    EdgeResMat[i][i] = GetResonanceIntegral(m.getBond(i));
                    if (EdgeResMat[i][i] == 0)
                        EdgeResMat[i][i] = 1;

                }

                DataMatrix = new Matrix(EdgeResMat);
            }
            
            // Calculates eigenvalues
            double[] eigenvalues;
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();

            Arrays.sort(eigenvalues);

            // SpMax
            this.SetByName("SpMax" + WEIGHT_SYMBOL[curWeight], eigenvalues[eigenvalues.length - 1]);
            
            // EEig
            for (int i=MinEig; i<=MaxEig; i++) {
                int idx = (eigenvalues.length - 1) - (i-1);
                if (idx>=0)
                    this.SetByName("EEig" + i + WEIGHT_SYMBOL[curWeight], eigenvalues[idx]);
                else
                    this.SetByName("EEig" + i + WEIGHT_SYMBOL[curWeight], 0);
            }

            // Spectral moment
            for (int i=MinSM; i<=MaxSM; i++) {
                double curSM = 0;
                double curSM_D5 = 0;
                for (int k=(eigenvalues.length-1); k>=0; k--) {
                    curSM += Math.pow(eigenvalues[k], (i));
                    if (Math.abs(eigenvalues[k]) > 0.000001)
                        curSM_D5 += Math.pow((eigenvalues[k]), (i));                    
                }
                curSM = Math.log(1 + curSM);
                curSM_D5 = Math.log(1 + curSM_D5);

                this.SetByName("ESpm" + i + WEIGHT_SYMBOL[curWeight], curSM);
                this.SetByName("ESpm" + i + WEIGHT_SYMBOL[curWeight] + "D5", curSM_D5);
            }      

        }        

    }

    
    
    private double GetResonanceIntegral(IBond bnd) {
        
        IAtom atA = bnd.getAtom(0);
        IAtom atB =  bnd.getAtom(1);
        String A = atA.getSymbol();
        String B = atB.getSymbol();
        
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("C") == 0)) )
            return 1.00;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("B") == 0)) ||
             ((A.compareToIgnoreCase("B") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("N") == 0)) ||
             ((A.compareToIgnoreCase("N") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.9;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("O") == 0)) ||
             ((A.compareToIgnoreCase("O") == 0) && (B.compareToIgnoreCase("C") == 0))) {
            if (bnd.getOrder() == IBond.Order.SINGLE)
                return 0.8;
            if (bnd.getOrder() == IBond.Order.DOUBLE)
                return 1.2;
        }
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("S") == 0)) ||
             ((A.compareToIgnoreCase("S") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("F") == 0)) ||
             ((A.compareToIgnoreCase("F") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.7;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Cl") == 0)) ||
             ((A.compareToIgnoreCase("Cl") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.4;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("Br") == 0)) ||
             ((A.compareToIgnoreCase("Br") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.3;
        if ( ((A.compareToIgnoreCase("C") == 0) && (B.compareToIgnoreCase("I") == 0)) ||
             ((A.compareToIgnoreCase("I") == 0) && (B.compareToIgnoreCase("C") == 0)))
            return 0.1;
        
        return 0.00;
    }
    
    
    private double GetDipoleMoment(IAtomContainer CurMol, IAtom at1, IAtom at2) {
        
        String a = at1.getSymbol();
        String b = at2.getSymbol();
        
        // C - something
        if (a.equalsIgnoreCase("C")) {
            
            // C-F
            if (b.equalsIgnoreCase("F")) {
                return 1.51;
            } 
            
            // C-Cl , C(Cl)-Cl , C(Cl)(Cl)-Cl
            if (b.equalsIgnoreCase("Cl")) {
                int nCl=0;
                for (IAtom at : CurMol.getConnectedAtomsList(at1)) {
                    if (at.getSymbol().equalsIgnoreCase("Cl"))
                        nCl++;
                }
                if (nCl==1)
                    return 1.56;
                if (nCl==2)
                    return 1.20;
                if (nCl==3)
                    return 0.83;
            } 
            
            // C-Br
            if (b.equalsIgnoreCase("Br")) {
                return 1.48;
            } 
        
            // C-I
            if (b.equalsIgnoreCase("I")) {
                return 1.29;
            } 
        
            // C-N , C=N , C#N
            if (b.equalsIgnoreCase("N")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.4;
                if (ord == IBond.Order.DOUBLE)
                    return 0.9;
                if (ord == IBond.Order.TRIPLE)
                    return 3.6;
            } 
        
            // C-O , C=O
            if (b.equalsIgnoreCase("O")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.86;
                if (ord == IBond.Order.DOUBLE)
                    return 2.4;
            } 
        
            // C-S , C=S
            if (b.equalsIgnoreCase("S")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 2.95;
                if (ord == IBond.Order.DOUBLE)
                    return 2.8;
            } 
            
        }
        

        // N-O , N-[O-] , N=O
        if ((a.equalsIgnoreCase("N")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nH=0;
            try { nH = at2.getImplicitHydrogenCount(); } catch (Exception e) {}
            int nConn = CurMol.getConnectedBondsCount(at2) + nH;
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 3.2;
            if ((ord == IBond.Order.SINGLE) && (nConn==1))
                return 2.0;
//                return 0.3;
            if ((ord == IBond.Order.DOUBLE) && (nConn==1))
                return 2.0;
        } 

        
        // S-[O-]
        if ((a.equalsIgnoreCase("S")) && (b.equalsIgnoreCase("O"))) {
            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nConn = CurMol.getConnectedBondsCount(at2);
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 2.9;
        } 

        
        // C(*)(*)-C(*)(*)(*) , C(*)(*)-C , CC(*)(*)(*)
        if ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("C"))) {

            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0;

            int nH1=0, nH2=0;
            try {
                nH1 = at1.getImplicitHydrogenCount();
            } catch (Exception E) {}
            try {
                nH2 = at2.getImplicitHydrogenCount();
            } catch (Exception E) {}
            
            int nConn1 = CurMol.getConnectedBondsCount(at1) + nH1;
            int nConn2 = CurMol.getConnectedBondsCount(at2) + nH2;

            if ((nConn1==3) && (nConn2==4))
                return 0.68;
            if ((nConn1==3) && (nConn2==2))
                return 1.15;
            if ((nConn1==2) && (nConn2==4))
                return 1.48;
        } 
        
        return 0;
    }
        

    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException 
     */
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        EdgeAdjacency block = new EdgeAdjacency();
        block.CloneDetailsFrom(this);
        return block;
    }

    
}
