package insilico.fish_nic_tk.descriptors.weights;

import Jama.Matrix;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;

import java.util.ArrayList;
import java.util.List;

/**
 * Walk And Path molecular descriptors.
 * Calculates MWC, MPC, PW, SRW
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Slf4j
public class WalkAndPath extends DescriptorBlock {

    private final static long serialVersionUID = 1L;
    private final static String BlockName = "Walk and Path Descriptors";

    public final static String PARAMETER_PATH_01 = "path01";
    public final static String PARAMETER_PATH_02 = "path02";
    public final static String PARAMETER_PATH_03 = "path03";
    public final static String PARAMETER_PATH_04 = "path04";
    public final static String PARAMETER_PATH_05 = "path05";
    public final static String PARAMETER_PATH_06 = "path06";
    public final static String PARAMETER_PATH_07 = "path07";
    public final static String PARAMETER_PATH_08 = "path08";
    public final static String PARAMETER_PATH_09 = "path09";
    public final static String PARAMETER_PATH_10 = "path10";
    
    public final static String PARAMETER_INCLUDE_SRW = "calcSRW";
    
    public final static String PARAMETER_INCLUDE_INDICES = "calcInd";

    private boolean defaultDescriptors;



    /**
     * Constructor. This should not be used, no weight is specified. The 
     * overloaded constructors should be used instead.
     */
    public WalkAndPath() {
        super();
        this.Name = WalkAndPath.BlockName;
        this.defaultDescriptors = true;
    }

    public WalkAndPath(boolean defaultDescriptors) {
        super();
        this.Name = WalkAndPath.BlockName;
        this.defaultDescriptors = defaultDescriptors;
    }

    
    
    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        ArrayList<Integer> pathsList = BuildPathsList();
        for (Integer curPath : pathsList) {
            Add("MWC" + curPath.toString(), "");
            Add("MPC" + curPath.toString(), "");
            Add("PW" + curPath.toString(), "");
            Add("piPC" + curPath.toString(), "");
        }
        if (getBoolProperty(PARAMETER_INCLUDE_SRW)) {
            int SRWMaxPath = pathsList.get(pathsList.size()-1);
            for (int i=1; i<=SRWMaxPath; i++)
                Add("SRW" + i, "");
        }
        if (getBoolProperty(PARAMETER_INCLUDE_INDICES)) {
            // if these indices are included, all paths from 1 to 9 are
            // automatically included and calculated
            Add("piID", "");
            Add("TPC", "");
            Add("PCR", "");
        }
        SetAllValues(Descriptor.MISSING_VALUE);
    }


    private ArrayList<Integer> BuildPathsList() {
        ArrayList<Integer> p = new ArrayList<>();
        if (defaultDescriptors) {
            setBoolProperty(PARAMETER_INCLUDE_INDICES, true);
            setBoolProperty(PARAMETER_INCLUDE_SRW, true);
        }
        boolean All = getBoolProperty(PARAMETER_INCLUDE_INDICES);
        if ((getBoolProperty(PARAMETER_PATH_01))||(All)) p.add(1);
        if ((getBoolProperty(PARAMETER_PATH_02))||(All)) p.add(2);
        if ((getBoolProperty(PARAMETER_PATH_03))||(All)) p.add(3);
        if ((getBoolProperty(PARAMETER_PATH_04))||(All)) p.add(4);
        if ((getBoolProperty(PARAMETER_PATH_05))||(All)) p.add(5);
        if ((getBoolProperty(PARAMETER_PATH_06))||(All)) p.add(6);
        if ((getBoolProperty(PARAMETER_PATH_07))||(All)) p.add(7);
        if ((getBoolProperty(PARAMETER_PATH_08))||(All)) p.add(8);
        if ((getBoolProperty(PARAMETER_PATH_09))||(All)) p.add(9);
        if ((getBoolProperty(PARAMETER_PATH_10))||(All)) p.add(10);

        return p;
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

        ArrayList<Integer> pathsList = BuildPathsList();
        if (pathsList.isEmpty())
            return;
        
        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] AdjMat = null;
        try {
            AdjMat = mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        int nSK = m.getAtomCount();
        double piID=0, TPC=0;
                
        // Cycle for all found path schemes
        for (Integer curPath : pathsList) {
            
            double CurPath=0, CurWalk=0, CurPW=0, CurMPC=0;
            int[] AtomWalk = GetAtomsWalks(curPath, AdjMatDbl);
            double[][] AtomPath = GetAtomsPaths(curPath, m);
            
            for (int i=0; i<nSK; i++) {
                CurWalk += AtomWalk[i];
                CurPath += AtomPath[i][0];
                CurMPC += AtomPath[i][1];
                CurPW += (double)AtomPath[i][0] / (double)AtomWalk[i];
            }
            
            if (curPath == 1) {
                CurWalk /= 2;
            } else {
                CurWalk = Math.log(1 + CurWalk);
            }
            CurPath /= 2;
            TPC += CurPath;
            CurMPC /= 2.0;
            piID += CurMPC;
            CurMPC = Math.log(1+CurMPC);
            CurPW /= nSK;
            
            
            SetByName("MWC" + curPath.toString(), CurWalk);
            SetByName("piPC" + curPath.toString(), CurMPC);
            SetByName("MPC" + curPath.toString(), CurPath);
            SetByName("PW" + curPath.toString(), CurPW);
            
        }
        
        if (getBoolProperty(PARAMETER_INCLUDE_INDICES)) {

            piID = Math.log(1+piID+nSK);
            TPC = Math.log(1+TPC+nSK);
            
            double PCR = piID / TPC;
            SetByName("piID", piID);
            SetByName("TPC", TPC);
            SetByName("PCR", PCR);
        }        
        
        if (getBoolProperty(PARAMETER_INCLUDE_SRW)) {
            int SRWMaxPath = pathsList.get(pathsList.size()-1);
            int[] MolSRW = GetSRWToLag(SRWMaxPath, AdjMatDbl);
            for (int i=1; i<=SRWMaxPath; i++)
                SetByName("SRW" + i, MolSRW[i]);
        }

    }

    
    
    private int[] GetAtomsWalks(int WalksOrder, double[][] AdjMatrix) {
        
        int[] walks = new int[AdjMatrix.length];
        
        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);
        
        for (int k=1; k<WalksOrder; k++) {
            mWalks = mWalks.times(mAdj);
        }
        
        for (int i=0; i<AdjMatrix.length; i++) {
            int CurSum = 0;
            for (int j=0; j<AdjMatrix[0].length; j++)
                CurSum += mWalks.get(i, j);
            walks[i] = CurSum;
        }
        
        return walks;
    }

    // index 0 = count of atom paths
    // index 1 = sum of paths weighted by bond order
    private double[][] GetAtomsPaths(int PathsOrder, IAtomContainer Mol) {
        
        int nSK = Mol.getAtomCount();
        double[][] paths = new double[nSK][2];
        
        for (int i=0; i<nSK; i++) {
            IAtom a = Mol.getAtom(i);
            int PathNum = 0;
            double totBO = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j)
                    continue;
                IAtom b =  Mol.getAtom(j);
                List<List<IAtom>> list=PathTools.getAllPaths(Mol, a, b);
                for (int k=0; k<list.size(); k++) {
                    if (list.get(k).size() == (PathsOrder+1)) {
                        PathNum++;
                        
                        double curBO=1;
                        for (int at_idx=0; at_idx<list.get(k).size()-1; at_idx++) {
                            IAtom at_1 = list.get(k).get(at_idx);
                            IAtom at_2 = list.get(k).get(at_idx+1);
                            IBond curBond = Mol.getBond(at_2, at_1);
                            if (curBond.getFlag(CDKConstants.ISAROMATIC)) 
                                curBO *= 1.5;
                            else {
                                if (curBond.getOrder() == Order.SINGLE) curBO *= 1;
                                if (curBond.getOrder() == Order.DOUBLE) curBO *= 2;
                                if (curBond.getOrder() == Order.TRIPLE) curBO *= 3;
                                if (curBond.getOrder() == Order.QUADRUPLE) curBO *= 4;
                            }
                        }
                        totBO+=curBO;
                    }
                }
                paths[i][0] = PathNum;
                paths[i][1] = totBO;
            }
        }
        
        return paths;
    }    
    
    
    private int[] GetSRWToLag(int Lag, double[][] AdjMatrix) {
        
        int nSK = AdjMatrix.length;
        int[] MolSRW = new int[Lag+1];
        
        Matrix mAdj = new Matrix(AdjMatrix);
        Matrix mWalks = new Matrix(AdjMatrix);
        
        for (int k=1; k<(Lag+1); k++) {
            
            MolSRW[k] = 0;
            for (int i=0; i<nSK; i++) {
                MolSRW[k] += mWalks.get(i, i);
            }
            
            mWalks = mWalks.times(mAdj);
        }
        
        return MolSRW;
    }
    
    

    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException 
     */
    @Override
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        WalkAndPath block = new WalkAndPath();
        block.CloneDetailsFrom(this);
        return block;
    }

    
}
