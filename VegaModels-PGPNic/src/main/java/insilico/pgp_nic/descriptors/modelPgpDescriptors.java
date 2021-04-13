
package insilico.pgp_nic.descriptors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.constant.InsilicoConstants;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.AtomCenteredFragments;
import insilico.core.descriptor.blocks.FunctionalGroups;
import insilico.core.descriptor.blocks.weights.TopologicalDistances;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;

import java.util.*;

import insilico.core.tools.utils.MoleculeUtilities;
import insilico.pgp_nic.descriptors.weights.*;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

/**
 *
 * @author Alberto
 */
public class modelPgpDescriptors {
    
    public double[] Descriptors;
    
    public modelPgpDescriptors() {
        Descriptors = new double[26];
        for (double d : Descriptors) 
            d = InsilicoConstants.MISSING_VALUE;
    }
    
    public void Calculate(InsilicoMolecule mol, FunctionalGroups desc_FC, 
                          AtomCenteredFragments desc_ACF) throws Exception {
        
        IAtomContainer curMol = mol.GetStructure();
        
        int nSK = curMol.getAtomCount();
        int nBO = curMol.getBondCount();
        int H[] = new int[nSK];
        int nTotH = 0;
        
        for (int i=0; i<nSK; i++) {

            IAtom CurAt =  curMol.getAtom(i);

            // Hydrogens
            H[i] = 0;
            try {
                H[i] = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) { }
            nTotH += H[i];
        }     
        
        // [0] H%
        Descriptors[0] = (nTotH/(double)(nSK + nTotH))*100;
        
        
        /////////////////
        
        
        int nR07 = 0;
        IRingSet allRings = mol.GetAllRings();
        Iterator<IAtomContainer> RingsIterator = allRings.atomContainers().iterator();        
        while (RingsIterator.hasNext()) {
        IRing ring = (IRing)RingsIterator.next();
            if (ring.getAtomCount() == 7)
                nR07++;
        }
        
        // [1] nR07
        Descriptors[1] = nR07;
        
        
        /////////////////
           
        
        double[][] DDrMatrix = mol.GetMatrixDistanceDetour();
        RingSet MolRings = mol.GetAllRings();

        double Ddtr11 = 0;

        for (int i=0; i< MolRings.getAtomContainerCount(); i++) {
            IRing r = (IRing) MolRings.getAtomContainer(i);
            int rSize = r.getAtomCount();

            if (rSize == 11) {
                for (IAtom at : r.atoms()) {
                    int atNum = curMol.getAtomNumber(at);
                    double rowSum = 0;
                    for (int k=0; k<DDrMatrix[atNum].length; k++)
                        rowSum += DDrMatrix[atNum][k];
                    Ddtr11  += rowSum;    
                }
            }

        }

        // [2] D/Dtr11
        Descriptors[2] = Ddtr11;
        
        
        /////////////////
        
        
        int[][] AdjMat = mol.GetMatrixAdjacency();
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];

        double CurWalk=0;
        int[] AtomWalk = GetAtomsWalks(1, AdjMatDbl);
        for (int i=0; i<nSK; i++) {
            CurWalk += AtomWalk[i];
        }
        CurWalk /= 2; // if (curPath == 1) 

        // [3] MWC01
        Descriptors[3] = Math.log(1 + CurWalk);
        
        
        /////////////////
        
        
        // Gets matrices
        double[][] ConnAugMatrix = mol.GetMatrixConnectionAugmented();
                
        int[] VD = VertexDegree.getWeights(curMol, true);
        double X2A = 0;
        
        // clears VD matrix from linked F
        for (int i=0; i<nSK; i++) 
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ((ConnAugMatrix[i][j]>0) && (ConnAugMatrix[j][j]==9))
                    VD[i]--;
            }
        
        int nPaths = 0;
        for (int i=0; i<nSK; i++) {

            if (ConnAugMatrix[i][i] == 9) 
                continue; // F not taken into account
                        
            IAtom at =  curMol.getAtom(i);
            List<List<IAtom>> CurPaths =  PathTools.getPathsOfLength(curMol, at, 2);
            nPaths += CurPaths.size();
            for (List<IAtom> curPath : CurPaths) {
                double prodX = 1;
                List<IAtom> CurPath = curPath;
                for (IAtom iAtom : CurPath) {
                    int atIdx = curMol.getAtomNumber(iAtom);
                    prodX *= VD[atIdx];
                }
                X2A += Math.pow(prodX, -0.5);
            }
        }
        
        X2A /= nPaths;
        
        // [4] X2A
        Descriptors[4] = X2A;
        
        
        /////////////////
        
        
        IAtomContainer mol_H = Manipulator.AddHydrogens(curMol);

        double[][] ConnMat_H = ConnectionAugMatrix.getMatrix(mol_H);
        int[][] TopoDistMat_H = TopoDistanceMatrix.getMatrix(mol_H);
        int nSK_H = mol_H.getAtomCount();
        int[] VertexDeg = VertexDegree.getWeights(mol_H, true);        
       
        ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK_H);
        for (int i=0; i<nSK_H; i++) {
            IAtom atStart = mol_H.getAtom(i);
            ArrayList<String> CurNeig = new ArrayList<>();
            for (int j=0; j<nSK_H; j++) {
                if (i==j) continue;
                if (TopoDistMat_H[i][j] == 3) {
                    IAtom atEnd = mol_H.getAtom(j);
//                    List<IAtom> sp = PathTools.getShortestPath(mol_H, atStart, atEnd);
                    ShortestPaths shortestPaths = new ShortestPaths(mol_H, atStart);
                    List<IAtom> sp = Arrays.asList(shortestPaths.atomsTo(atEnd));
                    StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                    for (int k=0; k<(sp.size()-1); k++) {
                        int a = mol_H.indexOf(sp.get(k));
                        int b = mol_H.indexOf(sp.get(k + 1));
                        if (ConnMat_H[a][b] == 1)
                            bufPath.append("s");
                        if (ConnMat_H[a][b] == 2)
                            bufPath.append("d");
                        if (ConnMat_H[a][b] == 3)
                            bufPath.append("t");
                        if (ConnMat_H[a][b] == 1.5)
                            bufPath.append("a");
                        bufPath.append(sp.get(k + 1).getAtomicNumber());
                        bufPath.append("(").append(VertexDeg[mol_H.getAtomNumber(sp.get(k + 1))]).append(")");
                    }
                    CurNeig.add(bufPath.toString());
                }
            }
            Collections.sort(CurNeig);
            NeigList.add(CurNeig);
        }        

        ArrayList<ArrayList<String>> G = new ArrayList<>();
        ArrayList<Integer> Gn = new ArrayList<>();
        for (int i=0; i<nSK_H; i++) {
            ArrayList<String> CurNeig = NeigList.get(i);
            boolean foundMatch = false;
            for (int k=0; k<G.size(); k++) {
                if (CompareNeigVector(CurNeig, G.get(k))) {
                    foundMatch = true;
                    int buf = Gn.get(k);
                    Gn.set(k, (buf+1));
                    break;
                }
            }
            if (!foundMatch) {
                G.add(CurNeig);
                Gn.add(1);
            }
        }        

        double ic=0;
        for (int i=0; i<Gn.size(); i++)
            ic += ((double)Gn.get(i)/nSK_H) * (Math.log((double)Gn.get(i)/nSK_H));
        ic = (-1.00 / Math.log(2)) * ic;

        // [5] SIC3
        // PICCOLE DISCREPANZE NEI VALORI
        Descriptors[5] = ic / Log(2, nSK_H);
        
        
        /////////////////
        
        
        double[][] Burden_S = new double[nSK][nSK];
        double[][] ConnMat = ConnectionAugMatrix.getMatrix(curMol);
        EState est = new EState(curMol);
        for (int i=0; i<nSK; i++) 
            for (int j=0; j<nSK; j++) {                
                if (i==j) {
                    Burden_S[i][j] = est.getIS()[i];
                } else {
                    if (ConnMat[i][j] > 0) {
                        Burden_S[i][j] = Math.sqrt(ConnMat[i][j]);
                        if ( (curMol.getConnectedAtomsCount(curMol.getAtom(i)) == 1) ||
                             (curMol.getConnectedAtomsCount(curMol.getAtom(j)) == 1))
                            Burden_S[i][j] += 0.1;
                    } else {
                        Burden_S[i][j] = 0.001;
                    }
                }                
            }        
        
        // Calculates eigenvalues
        Matrix DataMatrix = new Matrix(Burden_S);
        EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
        Matrix EigVectors = ed.getV();        
        
        double LastEigCoef = 0;
        for (int i=0; i<nSK; i++) {
            LastEigCoef += EigVectors.get(i, 0);
        }
        
        // [6] Ve1sign_B(s) - Coefficient sum of the last eigenvector from Burden matrix weighted by I-state
        Descriptors[6] = Math.abs(LastEigCoef);
        
        
        /////////////////
        
        
        double[] w_m = Mass.getWeights(mol_H);
        double[] w_v = VanDerWaals.getWeights(mol_H);
        EStateCorrectForH ES = new EStateCorrectForH(mol_H);
        double[] w_s = ES.getIS();
        for (int i=0; i<nSK_H; i++) {
            if (mol_H.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                w_s[i] = 1;
        }
             
        double wA_m = 0;        
        for (int i=0; i<nSK_H; i++) { if(w_m[i] == Descriptor.MISSING_VALUE) throw new Exception("Missing value in weights"); wA_m += w_m[i]; }
        wA_m = wA_m / ((double) nSK_H);
        double wA_v = 0;        
        for (int i=0; i<nSK_H; i++) { if(w_v[i] == Descriptor.MISSING_VALUE) throw new Exception("Missing value in weights"); wA_v += w_v[i]; }
        wA_v = wA_v / ((double) nSK_H);
        double wA_s = 0;        
        for (int i=0; i<nSK_H; i++) { if(w_s[i] == Descriptor.MISSING_VALUE) throw new Exception("Missing value in weights"); wA_s += w_s[i]; }
        wA_s = wA_s / ((double) nSK_H);
        
        // ATSC7m MATS6v GATS4s
        double ACS=0, MoranAC=0, GearyAC=0;       
        double denom_moran = 0, delta_moran = 0;        
        double denom_geary = 0, delta_geary = 0;        

        for (int i=0; i<nSK_H; i++) {

            denom_moran += Math.pow((w_v[i] - wA_v), 2);
            denom_geary += Math.pow((w_s[i] - wA_s), 2);

            for (int j=(i+1); j<nSK_H; j++) {
                if (TopoDistMat_H[i][j] == 7) {
                    ACS += Math.abs((w_m[i]-wA_m) * (w_m[j]-wA_m));
                }
                if (TopoDistMat_H[i][j] == 6) {
                    MoranAC += (w_v[i] - wA_v) * (w_v[j] - wA_v);
                    delta_moran++;
                }
                if (TopoDistMat_H[i][j] == 4) {
                    GearyAC += Math.pow((w_s[i] - w_s[j]), 2);
                    delta_geary++;
                }
            }
        }
        
        if (delta_moran > 0) {
            if (denom_moran == 0)
                MoranAC = 1;
            else
                MoranAC = (MoranAC / delta_moran) / (denom_moran / ((double)nSK_H));
        }
        
        if (delta_geary > 0) {
            if (denom_geary == 0)
                GearyAC = 0;
            else
                GearyAC = (GearyAC / (2.0 * delta_geary)) / (denom_geary / ((double)(nSK_H - 1)));
        }
        
        // PICCOLE DISCREPANZE NEI VALORI
        
        // [7] ATSC7m
        Descriptors[7] = ACS;

        // [8] MATS6v
        Descriptors[8] = MoranAC;

        // [9] GATS4s
        Descriptors[9] = GearyAC;
        
        
        /////////////////
        
        
        // Calculate VSA
        double[] VSA = new double[nSK_H];
        
        for (int i=0; i<nSK_H; i++) {

            IAtom at = mol_H.getAtom(i);

            double vdwR = GetVdWRadius(mol_H, at);
            if (vdwR == Descriptor.MISSING_VALUE) { VSA[i] = 0; continue; };

            double coef = 0;
            for (IAtom connAt : mol_H.getConnectedAtomsList(at)) {

                double connAt_vdwR = GetVdWRadius(mol_H, connAt);
                if (connAt_vdwR == Descriptor.MISSING_VALUE) { coef = Descriptor.MISSING_VALUE; break; };

                double refR = this.GetRefBondLength(mol_H, at, connAt);
                if (refR == Descriptor.MISSING_VALUE) { coef = Descriptor.MISSING_VALUE; break; };

                double g_1 = Math.max(  Math.abs(vdwR - connAt_vdwR) , refR ) ;
                double g_2 = vdwR + connAt_vdwR;
                double g  = Math.min(g_1, g_2) ;

                coef += ( Math.pow(connAt_vdwR,2) - Math.pow( (vdwR - g), 2) ) / (g) ;
            }

            if (coef == Descriptor.MISSING_VALUE)
                VSA[i] = 0;
            else
                VSA[i] = ( 4.0 * Math.PI * Math.pow(vdwR,2) ) - Math.PI * vdwR * coef ;
            if (VSA[i] < 0)
                VSA[i] = 0;            
        }        

        boolean typeD[] = new boolean[nSK_H];
        for (int i=0; i<nSK_H; i++) {
            boolean tD=false;
            IAtom at = mol_H.getAtom(i);
            int countH = 0;
            for (IAtom c : mol_H.getConnectedAtomsList(at))
                if (c.getSymbol().equalsIgnoreCase("H"))
                    countH++;
            
            if (at.getSymbol().equalsIgnoreCase("O")) 
                if ( (at.getFormalCharge() == 0) && (countH == 1))
                    tD = true;
            
            if (at.getSymbol().equalsIgnoreCase("N")) 
                if  ( (at.getFormalCharge() == 0) &&( (countH == 1) || (countH ==2) ) )
                    tD = true;
            
            typeD[i] = tD;
        }        
        
        double[] w = GhoseCrippenWeightsFixed.GetHydrophobiticty(mol_H);

        double PVSA_logP_3 = 0;
        double P_VSA_ppp_D = 0;
        for (int i=0; i<nSK_H; i++) {
            if (w[i] != Descriptor.MISSING_VALUE) 
                if ((w[i] >= -0.5) && (w[i] < -0.25) ) 
                    PVSA_logP_3 += VSA[i];
            if (typeD[i])
                P_VSA_ppp_D += VSA[i];
        }

        // [10] P_VSA_LogP_3
        Descriptors[10] = PVSA_logP_3;

        // [11] P_VSA_ppp_D
        Descriptors[11] = P_VSA_ppp_D;
                
        
        /////////////////
        
        DescriptorBlock FC;
        if (desc_FC == null) {
            FC = new FunctionalGroups();
            FC.Calculate(mol);
        } else
            FC = desc_FC;
        
        // [12] nRCOOR
        Descriptors[12] = FC.GetByName("nRCOOR").getValue();
        
        // [13] nArCONHR
        Descriptors[13] = FC.GetByName("nArCONHR").getValue();
        
        // [14] nArCO
        Descriptors[14] = FC.GetByName("nArCO").getValue();

        
        /////////////////
        
        
        DescriptorBlock ACF;
        if (desc_ACF == null) {
            ACF = new AtomCenteredFragments();
            ACF.Calculate(mol);
        } else
            ACF = desc_ACF;
        
        // [15] H-048
        Descriptors[15] = ACF.GetByName("H-048").getValue();
        
        
        /////////////////


        double SdsCH = 0;
        EState es = new EState(curMol);
        for (int at=0; at<curMol.getAtomCount(); at++) {
            
            IAtom curAt = curMol.getAtom(at);
            
            // Count H
            int nH = 0;
            try {
                nH = curAt.getImplicitHydrogenCount();
            } catch (Exception e) {
                nH = 0;
            }
            
            // formal charge
            int Charge;
            try {
                Charge = curAt.getFormalCharge();
            } catch (Exception e) {
                Charge = 0;
            }                    
            
            // Count bonds
            int nBnd=0, nSng = 0, nDbl = 0, nTri = 0, nAr=0;
            for (IBond b : curMol.getConnectedBondsList(curAt)) {
                if (b.getFlag(CDKConstants.ISAROMATIC)) {
                    nAr++;
                    nBnd++;
                    continue;
                }
                if (b.getOrder() == IBond.Order.SINGLE) {
                    nSng++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.DOUBLE) {
                    nDbl++;
                    nBnd++;
                }
                if (b.getOrder() == IBond.Order.TRIPLE) {
                    nTri++;
                    nBnd++;
                }
            }
            
            if (curAt.getSymbol().equalsIgnoreCase("C"))
                if ((nBnd == 2) && (nDbl == 1) && (nSng == 1) && (nH == 1))
                    SdsCH += es.getEState()[at];            
        }
                
        // [16] SdsCH
        Descriptors[16] = SdsCH;

        
        /////////////////

        
        int[][] TopoMat = mol.GetMatrixTopologicalDistance();
        
        ArrayList<String>[] AtomTypes = new ArrayList[nSK];
        for (int i=0; i<nSK; i++) {
        
            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt = curMol.getAtom(i);
            boolean tN=false, tP=false, tD=false, tL=false;
            
            // Hydrogens
            int mH = CurAt.getImplicitHydrogenCount();
            
            // [+]
            if (CurAt.getFormalCharge() > 0) {
                boolean NpOm = false;
                if (ConnAugMatrix[i][i] == 7) { 
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrix[i][j]==1) {
                            if (ConnAugMatrix[j][j] == 8) {
                                IAtom Oxy =  curMol.getAtom(j);
                                if (Oxy.getFormalCharge()!=0)
                                    NpOm = true;
                            }
                        }
                    }
                }
                if (!NpOm)
                    tP = true;
            }

            // [-]
            if (CurAt.getFormalCharge() < 0)
                tN = true;
            
            // O
            if (CurAt.getSymbol().equalsIgnoreCase("O")) {                
                if ( (CurAt.getFormalCharge() == 0) && (mH == 1))
                    tD = true;                
            }
                
            // N (NH2 and N without H)
            if (CurAt.getSymbol().equalsIgnoreCase("N")) {

                int nSglBnd = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrix[i][j]>0) {
                        if (ConnAugMatrix[i][j] == 1)
                            nSglBnd++;
                        else
                            nOtherBnd++;
                    }
                }
                
                if ( (CurAt.getFormalCharge() == 0) &&
                     (mH == 2) &&
                     (nSglBnd == 1) &&
                     (nOtherBnd == 0) )   
                    tP = true;
                
                if  ( (CurAt.getFormalCharge() == 0) &&( (mH == 1) || (mH ==2) ) )
                    tD = true;
                
            }
            
            // COOH, POOH, SOOH
            if ( ( (CurAt.getSymbol().equalsIgnoreCase("C")) ||
                   (CurAt.getSymbol().equalsIgnoreCase("S")) ||
                   (CurAt.getSymbol().equalsIgnoreCase("P")) ) &&
                  (CurAt.getFormalCharge() == 0) )  {
                
                int nSglBnd = 0, nDblO = 0, nSglOH = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrix[i][j]>0) {
                        if (ConnAugMatrix[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrix[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrix[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrix[i][j] == 2) && (ConnAugMatrix[j][j] == 8) )
                                nDblO++;
                            else
                                nOtherBnd++;
                        }
                    }
                }
                
                if ( (nSglBnd == 2) && (nSglOH == 1) && (nDblO == 1) && (nOtherBnd == 0) )
                    tN = true;
            }
            
            if (CurAt.getSymbol().equalsIgnoreCase("Cl")) 
                tL = true;
            
            if (CurAt.getSymbol().equalsIgnoreCase("Br")) 
                tL = true;

            if (CurAt.getSymbol().equalsIgnoreCase("I")) 
                tL = true;
            
            if (CurAt.getSymbol().equalsIgnoreCase("C")) {
                boolean connOnlyToSingleC = true;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] > 1.5) ) {
                            connOnlyToSingleC = false;
                            break;
                        }
                    }
                }                
                if (connOnlyToSingleC)
                    tL = true;
            }
            
            if (CurAt.getSymbol().equalsIgnoreCase("S")) {
                boolean connOnlyToSingleC = true;
                int nSingleC = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrix[i][j]>0) {
                        if ( (ConnAugMatrix[j][j] != 6) || (ConnAugMatrix[i][j] != 1) ) {
                            connOnlyToSingleC = false;
                            break;
                        } else {
                            nSingleC++;
                        }
                    }
                }                
                if ( (connOnlyToSingleC) && (nSingleC == 2) )
                    tL = true;
            }
            
            if (tN) AtomTypes[i].add("N");
            if (tP) AtomTypes[i].add("P");
            if (tD) AtomTypes[i].add("D");
            if (tL) AtomTypes[i].add("L");            
        }

        int CATS2D_01_DN = 0;
        int CATS2D_05_PP = 0;
        int CATS2D_02_PL = 0;

        for (int i=0; i<nSK; i++) {
            for (int j=0; j<nSK; j++) {
                if (i==j) continue;
                if ( (isIn(AtomTypes[i], "D")) && (isIn(AtomTypes[j], "N")) && (TopoMat[i][j] == 1) ) 
                    CATS2D_01_DN++;
                if ( (isIn(AtomTypes[i], "P")) && (isIn(AtomTypes[j], "P")) && (TopoMat[i][j] == 5) ) 
                    CATS2D_05_PP++;
                if ( (isIn(AtomTypes[i], "P")) && (isIn(AtomTypes[j], "L")) && (TopoMat[i][j] == 2) ) 
                    CATS2D_02_PL++;
            }
        }
        
        // CATS2D_01_DN
        Descriptors[17] = CATS2D_01_DN;

        // CATS2D_05_PP
        Descriptors[18] = CATS2D_05_PP;
        
        // CATS2D_02_PL
        Descriptors[19] = CATS2D_02_PL;
        
        
        /////////////////
        
        // 20-25 B e F
        
        TopologicalDistances TP = new TopologicalDistances();
        TP.Calculate(mol);
        
        // B07[O-F]
        Descriptors[20] = TP.GetByName("B7(O..F)").getValue();
        
        // F01[C-C]
        Descriptors[21] = TP.GetByName("F1(C..C)").getValue();
        
        // F02[C-O]
        Descriptors[22] = TP.GetByName("F2(C..O)").getValue();
        
        // F04[C-P]
        Descriptors[23] = TP.GetByName("F4(C..P)").getValue();
        
        // F04[C-Br]
        Descriptors[24] = TP.GetByName("F4(C..Br)").getValue();
        
        // F07[O-F]
        Descriptors[25] = TP.GetByName("F7(O..F)").getValue();
        
        
    }
    
    
    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
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
      
    private boolean CompareNeigVector(ArrayList<String> A, ArrayList<String> B) {
        
        // note: integer in the vectors should be already sorted 
        
        if (A.size()!=B.size())
            return false;
        else {
            for (int i=0; i<A.size(); i++)
                if (!A.get(i).equalsIgnoreCase(B.get(i)))
                    return false;
        }
        
        return true;
    }
    
    private boolean AtomCouple (IAtom at1, IAtom at2, String symbol1, String symbol2) {
        if ( (at1.getSymbol().equalsIgnoreCase(symbol1)) && (at2.getSymbol().equalsIgnoreCase(symbol2)))
            return true;
        if ( (at1.getSymbol().equalsIgnoreCase(symbol2)) && (at2.getSymbol().equalsIgnoreCase(symbol1)))
            return true;
        return false;
    }
        
    private double GetRefBondLength(IAtomContainer m, IAtom at1, IAtom at2) {
        double len = Descriptor.MISSING_VALUE;
        for (int i=0; i<RefBondLengths.length; i++) {
            if (AtomCouple(at1, at2, (String)RefBondLengths[i][0], (String)RefBondLengths[i][1])) {
                len = (Double)RefBondLengths[i][2];
                break;
            }
        }
        
        // c - correction
        double c = 0;
        if (len != Descriptor.MISSING_VALUE) {
            double bnd = MoleculeUtilities.Bond2Double(m.getBond(at1, at2));
            if (bnd == 1.5) c = 0.1;
            if (bnd == 2) c = 0.2;
            if (bnd == 3) c = 0.3;
        }
        
        return len - c;
    }    
    
    private double GetVdWRadius(IAtomContainer m, IAtom at) {
        
        String s = at.getSymbol();
        if (s.equalsIgnoreCase("C")) 
            return 1.950;
        if (s.equalsIgnoreCase("N")) 
            return 1.950;
        if (s.equalsIgnoreCase("F")) 
            return  1.496;
        if (s.equalsIgnoreCase("P")) 
            return  2.287;
        if (s.equalsIgnoreCase("S")) 
            return  2.185;
        if (s.equalsIgnoreCase("Cl")) 
            return  2.044;
        if (s.equalsIgnoreCase("Br")) 
            return  2.166;
        if (s.equalsIgnoreCase("I")) 
            return 2.358;
        
        if (s.equalsIgnoreCase("O")) {
            
            // oxide
            if (m.getConnectedAtomsList(at).size() == 1)
                return 1.810;
            
            // acid
            int H = 0;
            for (IAtom c : m.getConnectedAtomsList(at))
                if (c.getSymbol().equalsIgnoreCase("H"))
                    H++;
            if (H > 0)
                return 2.152;
            
            return 1.779;
        }
        
        if (s.equalsIgnoreCase("H")) {
            IAtom connAt = m.getConnectedAtomsList(at).get(0);
            if (connAt.getSymbol().equalsIgnoreCase("O"))
                return 0.8;
            if (connAt.getSymbol().equalsIgnoreCase("N"))
                return 0.7;
            if (connAt.getSymbol().equalsIgnoreCase("P"))
                return 0.7;
            return 1.485;
        }
        
        return Descriptor.MISSING_VALUE;
    }        
    
    private static double Log(int base,double x) {
        double Logbx = Math.log10(x)/Math.log10((double)base);
        return Logbx;
    }        
    
    private final static Object[][] RefBondLengths = {
        {"Br","Br",2.54},
        {"Br","C",1.97},
        {"Br","Cl",2.36},
        {"Br","F",1.85},
        {"Br","H",1.44},
        {"Br","I",2.65},
        {"Br","N",1.84},
        {"Br","O",1.58},
        {"Br","P",2.37},
        {"Br","S",2.21},
        {"C","C",1.54},
        {"C","Cl",1.8},
        {"C","F",1.35},
        {"C","H",1.06},
        {"C","I",2.12},
        {"C","N",1.47},
        {"C","O",1.43},
        {"C","P",1.85},
        {"C","S",1.81},
        {"Cl","Cl",2.31},
        {"Cl","F",1.63},
        {"Cl","H",1.22},
        {"Cl","I",2.56},
        {"Cl","N",1.74},
        {"Cl","O",1.41},
        {"Cl","P",2.01},
        {"Cl","S",2.07},
        {"F","F",1.28},
        {"F","H",0.87},
        {"F","I",2.04},
        {"F","N",1.41},
        {"F","O",1.32},
        {"F","P",1.5},
        {"F","S",1.64},
        {"H","I",1.63},
        {"H","N",1.01},
        {"H","O",0.97},
        {"H","P",1.41},
        {"H","S",1.31},
        {"I","I",2.92},
        {"I","N",2.26},
        {"I","O",2.14},
        {"I","P",2.49},
        {"I","S",2.69},
        {"N","N",1.45},
        {"N","O",1.46},
        {"N","P",1.6},
        {"N","S",1.76},
        {"O","O",1.47},
        {"O","P",1.57},
        {"O","S",1.57},
        {"P","P",2.26},
        {"P","S",2.07},
        {"S","S",2.05}        
    };    
}
