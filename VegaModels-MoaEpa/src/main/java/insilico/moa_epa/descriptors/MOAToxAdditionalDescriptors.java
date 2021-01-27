package insilico.moa_epa.descriptors;

import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.AtomicNumber;
import insilico.core.molecule.tools.Manipulator;
import insilico.moa_epa.descriptors.weights.Mass;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.Hashtable;
import java.util.List;

import static insilico.moa_epa.descriptors.weights.ValenceVertexDegree.GetValenceElectronsNumber;


/**
 *
 * @author Alberto
 */
public class MOAToxAdditionalDescriptors {
    
    private Hashtable<String, Double> Descriptors;
    private Pattern[] q;
    
    private final String[][] SMARTS = {
        { "SdsssP_acnt", "[P](-*)(-*)(-*)(=*)" },
        { "-O- [2 phosphorus attach]", "O([P])[P]" },
        { "-OH [phosphorus attach]", "[O;D1]-[P]" },
        { "-NH- [aliphatic attach]", "[N;D2]([C])[C]" },
        { "-NH- [aromatic attach]", "[N;D2]([c])" },
        { "-C(=O)O- [nitrogen attach]", "C([N])(=O)[O;D2]" },
        { "-CH< [aromatic attach]", "[C;D3](-*)(-*)[a]" },
        { ">C= [aromatic attach]", "[C](=C)(-*)[a]" },
        { "C=O(ketone, aliphatic attach)", "C(=O)([C;!$(C=*)])[C;!$(C=*)]" },
        { "-C(=O)- [aromatic attach]", "C(=O)([c])C" },
        { "-C(=O)O- [cyclic]", "[C;R](=O)([O;D2])-*" },
        { "-C#N [aliphatic nitrogen attach]", "C(#N)N" },
        { "-NO2 [nitrogen attach]", "[N+](=O)([O-])N" },
        { "SdCH2_acnt", "[C;D1;H2]=*"},
        { "StsC_acnt", "C(-*)#*"},
        { "SdsssP_acnt", "P(-*)(-*)(-*)=*"},
        { "-NH- [nitrogen attach]", "[N;D2;H]([N])-*"},
        { "-CHO [aliphatic attach]", "[C;D2;H](=O)C"},
        { "-C(=O)O- [nitrogen attach]", "C(N)(=O)[O;D2]"},
        { "-C(=O)O- [aliphatic attach]", "C([C;!$(C=*)])(=O)[O;D2]"},
        { "StCH_acnt", "[C;D1;H]#*"},
        { "-S- [sulfur attach]", "[S;D2](-*)-[S]"},
        { "-CHO [aromatic attach]", "[C;D2;H](=O)[c]"},
        { "CH2=CHC(=O)O-", "[C;D1;H2]=[C;D2;H]C(=O)[O;D2]"},
        { "-Br [aromatic attach]", "[Br]-[a]"},
        { "-Cl [aromatic attach]", "[Cl]-[a]"},
        { "-I [aromatic attach]", "[I]-[a]"},
        { "-OH [aromatic attach]", "[O;D1;H]-[a]"}           
    };
    

    
    
    public MOAToxAdditionalDescriptors() throws InitFailureException {

        Descriptors = new Hashtable<>();
        
        int n = SMARTS.length;
        q = new Pattern[n];
        for (int i=0; i<n; i++)
            try {
                    q[i] = SmartsPattern.create(SMARTS[i][1]).setPrepare(false);
            } catch (Exception ex) {
                throw new InitFailureException("unable to init SMARTS parsers with: " + SMARTS[i][1]);
            }
    }    

    
    public void Calculate(InsilicoMolecule mol) throws GenericFailureException {
        
        //// SMARTS
        
        int n = SMARTS.length;
        
        for (int i=0; i<n; i++)
            Descriptors.put(SMARTS[i][0], 0.0);
                
        try {

            int mappings = 0;
            boolean status = false;
//            CustomQueryMatcher Matcher = new CustomQueryMatcher(mol);
            
            for (int i=0; i<n; i++) {
                if (q[i].matches(mol.GetStructure())) {
                    mappings = q[i].matchAll(mol.GetStructure()).countUnique();
                    Descriptors.put(SMARTS[i][0], (double) mappings);
                }
            }
            
        } catch (Exception e) {
            throw new GenericFailureException("unable to perform matching - " + e.getMessage());
        }                
        
        
        //// AMW (not scaled)
        
        double AMW = 0;
        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
            int nSK = curMol.getAtomCount();
            int nTotH=0;
            int H[] = new int[nSK];
            for (int i=0; i<nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];
            }
            
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");
            double[] wMassNS = new double[wMass.length];
            double CarbonWeight = 12.011;
            double HMassNS = HMass * CarbonWeight;                
            for (int i=0; i<nSK; i++)
                wMassNS[i] = wMass[i] * CarbonWeight;

            for (int i=0; i<nSK; i++) {
                AMW += wMassNS[i];
                if (H[i]>0) 
                    AMW += HMassNS * H[i];
            }            
            
            AMW = AMW / (nSK + nTotH);
                    
        } catch (InvalidMoleculeException e) {
            AMW = 0;
        }
        
        Descriptors.put("AMW", AMW);
        

        //// icycem
        
        double icycem = 0;
        try {            
            double icyce;
            RingSet rings = mol.GetAllRings();
            int C = rings.getAtomContainerCount();
            int A = mol.GetStructure().getAtomCount();

            if (C>0) {

                int n0=0;		
                int n1=0;

                for (int i = 0; i < rings.getAtomContainerCount(); i++) {
                    int ringSizei=rings.getAtomContainer(i).getAtomCount();
                    n1+=ringSizei;
                    n0+=A-ringSizei;			
                }

                if (n0>0) 
                    icyce=A*C*Log(2,A*C)-n1*Log(2,n1)-n0*Log(2,n0);
                else 
                    icyce=A*C*Log(2,A*C)-n1*Log(2,n1);

                icycem=icyce/(A*C);

            } else {
                icycem=0;
            }
        
        } catch (InvalidMoleculeException e) {
            icycem = 0;
        }
    
        Descriptors.put("icycem", icycem);
        
        
        //// IC

        try {
            double IC = IC(mol);
            Descriptors.put("ic", IC);
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate IC - " + e.getMessage());
        }      
        
        
        //// Chi chain 10
        try {
            double Xch10 = ChiChain(mol, 10);
            Descriptors.put("Xch10", Xch10);
        } catch (Exception e) {
            throw new GenericFailureException("unable to calculate Xch10 - " + e.getMessage());
        }      
        
        
    }
    
    private static double Log(int base,double x) {
        double Logbx = Math.log10(x)/ Math.log10((double)base);
        return Logbx;
    }        
    
    
    public double GetByName(String Name) throws DescriptorNotFoundException {
        if (!Descriptors.containsKey(Name))
            throw new DescriptorNotFoundException();
        return Descriptors.get(Name);
    }

    
    private double IC(InsilicoMolecule mol) throws Exception {
        
        int nSK = mol.GetStructure().getAtomCount();
        double[][] ValenceTopoMatrix = new double[nSK][nSK];

        AtomicNumber ZFinder;
        ZFinder = new AtomicNumber();
                
        for (int i=0; i<(nSK-1); i++) {
            for (int j=i+1; j<nSK; j++) {
                IAtom at_i =  mol.GetStructure().getAtom(i);
                IAtom at_j =  mol.GetStructure().getAtom(j);
                List<IAtom> path = PathTools.getShortestPath(mol.GetStructure(), at_i, at_j);
                
                // geometric mean of valence vertex degree of atoms in path
                double GM = 1;
                int n = path.size(); // no. of atoms in path
                for (IAtom at : path) {
                    
                    // calculate Valence Vertex Degree
                    int Z = ZFinder.GetAtomicNumber(at.getSymbol());
                    int Zv = GetValenceElectronsNumber(at.getSymbol());
                    int h = Manipulator.CountImplicitHydrogens(at);
                    int ch = at.getFormalCharge();
                    if (Zv == -999)
                        throw new Exception("can't calculate Zv");
                    double ValenceVertexDegree = (double)(Zv - h - ch) / (double)(Z - Zv - 1.00);
                    
                    GM *= ValenceVertexDegree;
                }                
                GM = Math.pow(GM, (1.0 / (double) n));
                
                ValenceTopoMatrix[i][j] = GM * (1.0 / (double) n);
                ValenceTopoMatrix[j][i] = ValenceTopoMatrix[i][j];
            }
        }
        
        // calculate Sv for all atoms
        double[] Sv = new double[nSK];
        for (int i=0; i<nSK; i++) {
            Sv[i] = 0;
            for (int j=0; j<nSK; j++) {
                if (i==j) 
                    continue;
                Sv[i] += ValenceTopoMatrix[i][j];
            }
        }
        
        // find eq classes
        double tol = 1.0E-6D;
        
        Hashtable<Double, Integer> EqClasses = new Hashtable<>();
        for (int i = 0; i <nSK; i++) {
            boolean foundKey = false;
            for (Double CurSv : EqClasses.keySet()) {
                if (Math.abs(Sv[i] - CurSv) < tol) {
                    foundKey = true;
                    Integer counter = EqClasses.get(CurSv);
                    counter++;
                    EqClasses.put(CurSv, counter);
                }
            }
            if (!foundKey)
                EqClasses.put(Sv[i], 1);
        }
        
        // calculate IC
        double IC = 0.0;
        for (Double CurSv : EqClasses.keySet()) {
            Integer count = EqClasses.get(CurSv);
            IC += count * Log(2, count);    
        }
        
        return IC;
    }
    
    
    private double ChiChain(InsilicoMolecule mol, int Lag) throws Exception {
        
        double Chi = 0;
        RingSet rs = mol.GetAllRings();
        
        for (IAtomContainer curRing : rs.atomContainers()) {
            if (curRing.getAtomCount() == Lag) {
                double c = 1;
                for (IAtom a : curRing.atoms()) {
                    int vertexDegree = mol.GetStructure().getConnectedBondsCount(a);
                    c *= Math.pow(vertexDegree, -0.5);
                }
                Chi += c;
            }            
        }
                
        return Chi;
    }    
}
