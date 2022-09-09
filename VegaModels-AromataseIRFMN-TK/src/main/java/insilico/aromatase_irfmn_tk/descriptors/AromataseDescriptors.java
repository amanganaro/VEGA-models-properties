package insilico.aromatase_irfmn_tk.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.aromatase_irfmn_tk.descriptors.weights.Mass;
import insilico.core.descriptor.Descriptor;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.log4j.Log4j;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Arrays;
import java.util.List;

/**
 * Edge Adjacency descriptors (on augmented matrix).
 * 
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */

public class AromataseDescriptors {

    private static final Logger log = LogManager.getLogger(AromataseDescriptors.class);

    public static double Calculate_SM6_Bm(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return (Descriptor.MISSING_VALUE);
        }
        
        int nSK = m.getAtomCount();
        
        // Gets matrices
        double[][] BurdenMat;
        try {
            BurdenMat = mol.GetMatrixBurden();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return(Descriptor.MISSING_VALUE);
        }
        
        double[] w = Mass.getWeights(m);
        
        // If one or more weights are not available, sets all to missing value
        boolean MissingWeight = false;
        for (int i=0; i<nSK; i++) 
            if (w[i] == Descriptor.MISSING_VALUE)
                MissingWeight = true;
        if (MissingWeight)        
            return Descriptor.MISSING_VALUE;

        // Builds the weighted matrix
        for (int i=0; i<nSK; i++) {
            BurdenMat[i][i] = w[i];
        }

        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(BurdenMat);
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();
        Arrays.sort(eigenvalues);
        
        double sum = 0;
        for (double val : eigenvalues)
            sum += Math.pow(val, 6);

        return Math.log(1 + sum);
    }
    
    
    public static double Calculate_spMaxAEAdm(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return (Descriptor.MISSING_VALUE);
        }
        
        int nBO = m.getBondCount();
        
        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            return(Descriptor.MISSING_VALUE);
        }
        
        // Gets matrices
        double[][][] EdgeAdjMat = null;
        try {
            EdgeAdjMat = mol.GetMatrixEdgeAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            return(Descriptor.MISSING_VALUE);
        }


        // Weight dm 
        double[][] EdgeDipoleMat = new double[EdgeAdjMat.length][EdgeAdjMat[0].length];
        for (int i=0; i<EdgeAdjMat.length; i++) {
            for (int j=0; j<EdgeAdjMat[0].length; j++) {
                EdgeDipoleMat[i][j] = EdgeAdjMat[i][j][0];
            }
        }

        for (int i=0; i<m.getBondCount(); i++) {
            IAtom a =  m.getBond(i).getAtom(0);
            IAtom b =  m.getBond(i).getAtom(1);
            double CurVal = GetDipoleMoment(m, a, b);
            if (CurVal == 0)
                CurVal = GetDipoleMoment(m, b, a);
            EdgeDipoleMat[i][i] = CurVal;     
        }
        
        Matrix DataMatrix = new Matrix(EdgeDipoleMat);

        // Calculates eigenvalues
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        // SpMax
        double spMax = eigenvalues[eigenvalues.length - 1];

        return spMax / (double) nBO;
    }

        
    public static double Calculate_SM15_EAed(InsilicoMolecule mol) {

        IAtomContainer m;
        try {
            m = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            return (Descriptor.MISSING_VALUE);
        }
        
        int nBO = m.getBondCount();
        
        // Only for mol with nSK>1
        if (m.getAtomCount() < 2) {
            return(Descriptor.MISSING_VALUE);
        }
        
        // Gets matrices
        double[][] EdgeAdjMat = new double[nBO][nBO];
        for (int i=0; i<nBO; i++) {
            IBond b = m.getBond(i);
            for (int k=0; k<b.getAtomCount(); k++) {
                IAtom at = b.getAtom(k);
                List<IBond> connBonds = m.getConnectedBondsList(at);
                for (int z=0; z<m.getConnectedBondsCount(at); z++) {
                    IBond connBond = connBonds.get(z);
                    int BondNum = m.getBondNumber(connBond);
                    if (BondNum != i) {
                        EdgeAdjMat[i][BondNum] = 1;
                        EdgeAdjMat[BondNum][i] = 1;
                    }
                }           
            }
        }
        double[] wEdgeDegree = new double[nBO];
        for (int i=0; i<nBO; i++) {
            int ed = 0;
            for (int j=0; j<nBO; j++)
                ed += EdgeAdjMat[i][j];
            wEdgeDegree[i] = ed;
        }        
        double[][] EdgeAdjMatEd = new double[nBO][nBO];
        for (int i=0; i<nBO; i++) {
            for (int j=0; j<nBO; j++) {
                if (EdgeAdjMat[i][j] == 1) {
                    EdgeAdjMatEd[i][j] = wEdgeDegree[j];
                    EdgeAdjMatEd[j][i] = wEdgeDegree[i];
                } else {
                    EdgeAdjMatEd[i][j] = 0;
                    EdgeAdjMatEd[j][i] = 0;                    
                }
            }
        }
        
        Matrix DataMatrix = new Matrix(EdgeAdjMatEd);
        
        // Calculates eigenvalues
        double[] eigenvalues;
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        eigenvalues = ed.getRealEigenvalues();

        Arrays.sort(eigenvalues);

        double sum = 0;
        for (double val : eigenvalues)
            sum += Math.pow(val, 15);

        return Math.log(1 + sum);
    }

    
    private static double GetDipoleMoment(IAtomContainer CurMol, IAtom at1, IAtom at2) {
        
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
                    return 0.0;
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
                    return 0.0;
                IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
                if (ord == IBond.Order.SINGLE)
                    return 0.86;
                if (ord == IBond.Order.DOUBLE)
                    return 2.4;
            } 
        
            // C-S , C=S
            if (b.equalsIgnoreCase("S")) {
                if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                    return 0.0;
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
                return 0.0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nH=0;
            try { nH = at2.getImplicitHydrogenCount();
            } catch (Exception e) {
                log.warn(e.getMessage());
            }
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
                return 0.0;
            IBond.Order ord = CurMol.getBond(at1, at2).getOrder();
            int nConn = CurMol.getConnectedBondsCount(at2);
            if ((ord == IBond.Order.SINGLE) && (nConn==2))
                return 2.9;
        } 

        
        // C(*)(*)-C(*)(*)(*) , C(*)(*)-C , CC(*)(*)(*)
        if ((a.equalsIgnoreCase("C")) && (b.equalsIgnoreCase("C"))) {

            if (CurMol.getBond(at1, at2).getFlag(CDKConstants.ISAROMATIC))
                return 0.0;

            int nH1=0, nH2=0;
            try {
                nH1 = at1.getImplicitHydrogenCount();
            } catch (Exception E) {
                log.warn(E.getMessage());
            }
            try {
                nH2 = at2.getImplicitHydrogenCount();
            } catch (Exception E) {
                log.warn(E.getMessage());
            }
            
            int nConn1 = CurMol.getConnectedBondsCount(at1) + nH1;
            int nConn2 = CurMol.getConnectedBondsCount(at2) + nH2;

            if ((nConn1==3) && (nConn2==4))
                return 0.68;
            if ((nConn1==3) && (nConn2==2))
                return 1.15;
            if ((nConn1==2) && (nConn2==4))
                return 1.48;
        } 
        
        return 0.0;
    }
        
    
}
