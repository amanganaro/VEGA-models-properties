package insilico.skin_irfmn.desciptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.ConnectionAugMatrix;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.MoleculeUtilities;
import insilico.skin_irfmn.desciptors.weights.DescriptorMLogP;
import insilico.skin_irfmn.desciptors.weights.EStateCorrectForH;
import insilico.skin_irfmn.desciptors.weights.FG;
import insilico.skin_irfmn.desciptors.weights.VertexDegree;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ShortestPaths;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Slf4j
@Data
public class EmbeddedDescriptors {

    private final double MISSING_VALUE = -999;

    private double nDblBo = MISSING_VALUE;
    private double nRNO2 = MISSING_VALUE;
    private double nArNO2 = MISSING_VALUE;
    private double IC1 = MISSING_VALUE;
    private double SIC2 = MISSING_VALUE;
    private double GATS8s = MISSING_VALUE;
    private double nRCHO = MISSING_VALUE;
    private double CATS2D_1_NL = MISSING_VALUE;
    private double F4_C_O = MISSING_VALUE;
    private double MLOGP = MISSING_VALUE;
    
    private final static String TYPE_N = "N";
    private final static String TYPE_L = "L";

    private final static String[][] AtomCouplesCATS = {
            {TYPE_N, TYPE_L},
    };

    private double[][] ConnAugMatrixCATS;
    

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol){
        CalculateFG(Mol);
        CalculateMLogP(Mol);
        CalculateConstitutional(Mol);
        CalculateInformationCOntent(Mol);
        CalculateACH(Mol);
        CalculateCATS2D(Mol);
        CalculateTopologicalDistances(Mol);
    }

    private void CalculateTopologicalDistances(InsilicoMolecule Mol){

        this.setF4_C_O(0);
        String[][] AtomCouples = {
                {"C", "O"},
        };

        int MAX_TOPO_DISTANCE = 10;

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setF4_C_O(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setF4_C_O(MISSING_VALUE);
            return;
        }

        for (int d=0; d<AtomCouples.length; d++) {

            int[] descF = new int[MAX_TOPO_DISTANCE];
            Arrays.fill(descF, 0);

            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase(AtomCouples[d][0])) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if (m.getAtom(j).getSymbol().equalsIgnoreCase(AtomCouples[d][1])) {

                            if (TopoMat[i][j] <= MAX_TOPO_DISTANCE) {
                                descF[TopoMat[i][j]-1]++;
                            }

                        }
                    }
                }
            }

            // Fix: if atoms are the same, resulting value is calculated twice
            if (AtomCouples[d][0].compareTo(AtomCouples[d][1]) == 0) {
                for (int i=0; i<descF.length; i++)
                    descF[i] /= 2;
            }


            for (int i=0; i<descF.length; i++){
                if(i == 3){
                    this.setF4_C_O(descF[i]);
                    return;
                }

            }
        }


    }

    private void CalculateCATS2D(InsilicoMolecule Mol){
        this.setCATS2D_1_NL(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setCATS2D_1_NL(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Gets matrices
        int[][] TopoMat = null;
        try {
            TopoMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_1_NL(MISSING_VALUE);
            return;
        }
        try {
            ConnAugMatrixCATS = Mol.GetMatrixConnectionAugmented();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setCATS2D_1_NL(MISSING_VALUE);
            return;
        }

        // Gets CATS types
        ArrayList<String>[] CatsTypes = setCatsAtomType(m);


        for (int d = 0; d< AtomCouplesCATS.length; d++) {

            int descT = 0;
            int[] desc = new int[10];
            Arrays.fill(desc, 0);

            for (int i=0; i<nSK; i++) {
                if ( isIn(CatsTypes[i], AtomCouplesCATS[d][0]) ) {
                    for (int j=0; j<nSK; j++) {
                        if (i==j) continue;
                        if ( isIn(CatsTypes[j], AtomCouplesCATS[d][1]) ) {

                            if (TopoMat[i][j] < 10)
                                desc[TopoMat[i][j]]++;

                        }
                    }
                }
            }

            for (int i=0; i<desc.length; i++)
                if(i == 1){
                    this.setCATS2D_1_NL(desc[i]);
                    return;
                }
        }
        
    }

    private void CalculateACH(InsilicoMolecule Mol){
        this.setGATS8s(0);

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            this.setGATS8s(MISSING_VALUE);
            return;
        }

        // Gets matrices
        int[][] TopoMatrix;
        try {
            TopoMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            this.setGATS8s(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();
        double[] w;
        try {
            EStateCorrectForH ES = new EStateCorrectForH(m);
            w = ES.getIS();

            // correction for compatibility with D7
            // H I-state is always 1
            for (int i=0; i<nSK; i++) {
                if (m.getAtom(i).getSymbol().equalsIgnoreCase("H"))
                    w[i] = 1;
            }
        } catch (Exception e) {
            w = new double[nSK];
            for (int i=0; i<nSK; i++) w[i]=Descriptor.MISSING_VALUE;
        }

        double wA = 0;
        for (int i=0; i<nSK; i++)
            wA += w[i];
        wA = wA / ((double) nSK);

        double GearyAC=0;
        double denom = 0, delta = 0;

        for (int i=0; i<nSK; i++) {

            denom += Math.pow((w[i] - wA), 2);

            for (int j=0; j<nSK; j++)
                if (TopoMatrix[i][j] == 8) {
                    GearyAC += Math.pow((w[i] - w[j]), 2);
                    delta++;
                }
        }

        if (delta > 0) {
            if (denom == 0) {
                GearyAC = 0;
            } else {
                GearyAC = ((1 / (2 * delta)) * GearyAC) / ((1 / ((double)(nSK - 1))) * denom);
            }
        }

        this.setGATS8s(GearyAC);




    }

    private void CalculateInformationCOntent(InsilicoMolecule Mol){
        this.setIC1(0);
        this.setSIC2(0);

        int MaxPath = 2;

        IAtomContainer m;
        try {
            IAtomContainer orig_m = Mol.GetStructure();
            m = Manipulator.AddHydrogens(orig_m);
        } catch (InvalidMoleculeException | GenericFailureException e) {
            this.setIC1(MISSING_VALUE);
            this.setSIC2(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] ConnMat;
        int[][] TopoDistMat;
        try {
            ConnMat = ConnectionAugMatrix.getMatrix(m);
            TopoDistMat = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            this.setIC1(MISSING_VALUE);
            this.setSIC2(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);


        // Neighborhood indices

        double bic_denom = 0;
        for (IBond bnd : m.bonds())
            bic_denom += MoleculeUtilities.Bond2Double(bnd);
        bic_denom = Math.log(bic_denom) / Math.log(2);

        for (int CurLag=1; CurLag<=MaxPath; CurLag++) {

            // Create belonging class for each atom(vertex)
            ArrayList<ArrayList<String>> NeigList = new ArrayList<>(nSK);
            for (int i=0; i<nSK; i++) {
                IAtom atStart =  m.getAtom(i);
                ArrayList<String> CurNeig = new ArrayList<>();
                for (int j=0; j<nSK; j++) {
                    if (i==j) continue;
                    if (TopoDistMat[i][j] == CurLag) {
                        IAtom atEnd =  m.getAtom(j);
                        ShortestPaths shortestPaths = new ShortestPaths(m, atStart);
                        List<IAtom> sp = Arrays.asList(shortestPaths.atomsTo(atEnd));
                        // DEPRECATED METHOD
//                        List<IAtom> sp = PathTools.getShortestPath(m, atStart, atEnd);

                        StringBuilder bufPath = new StringBuilder("" + sp.get(0).getAtomicNumber());
                        for (int k=0; k<(sp.size()-1); k++) {
                            int a = m.indexOf(sp.get(k));
                            int b = m.indexOf(sp.get(k + 1));
                            if (ConnMat[a][b] == 1)
                                bufPath.append("s");
                            if (ConnMat[a][b] == 2)
                                bufPath.append("d");
                            if (ConnMat[a][b] == 3)
                                bufPath.append("t");
                            if (ConnMat[a][b] == 1.5)
                                bufPath.append("a");
                            bufPath.append(sp.get(k + 1).getAtomicNumber());
                            bufPath.append("(").append(VertexDeg[m.indexOf(sp.get(k + 1))]).append(")");
                        }
                        CurNeig.add(bufPath.toString());
                    }
                }
                Collections.sort(CurNeig);
                NeigList.add(CurNeig);
            }

            // Calculates equivalence classes
            ArrayList<ArrayList<String>> G = new ArrayList<>();
            ArrayList<Integer> Gn = new ArrayList<>();
            for (int i=0; i<nSK; i++) {
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

            // Calculate IC and CIC indices

            double ic=0;
            for (int i=0; i<Gn.size(); i++)
                ic += ((double)Gn.get(i)/nSK) * (Math.log((double)Gn.get(i)/nSK));
            ic = (-1.00 / Math.log(2)) * ic;

            double diff = Math.log(nSK) / Math.log(2);

            double cic = diff - ic;

            double bic = bic_denom==0 ? 0 : ic / bic_denom;

            double sic = ic / diff;

            if(CurLag == 1)
                this.setIC1(ic);

            if(CurLag == 2)
                this.setSIC2(sic);

        }
    }

    private void CalculateConstitutional(InsilicoMolecule Mol){
        this.setNDblBo(0);

        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setNDblBo(MISSING_VALUE);
            return;
        }

        try {

            int nBO = curMol.getBondCount();

            int nArBonds = 0, nDblBonds = 0, nTrpBonds = 0, nMulBonds = 0;
            double scbo = 0;

            for (int i = 0; i < nBO; i++) {

                IBond CurBo = curMol.getBond(i);

                if (CurBo.getFlag(CDKConstants.ISAROMATIC)) {
                    nArBonds++;
                    nMulBonds++;
                    scbo += 1.5;
                } else {
                    if (CurBo.getOrder() == IBond.Order.SINGLE) {
                        scbo++;
                    } else {
                        nMulBonds++;
                        if (CurBo.getOrder() == IBond.Order.DOUBLE) {
                            nDblBonds++;
                            scbo += 2;
                        }
                    }
                }

            }
            this.setNDblBo(nDblBonds);



        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setNDblBo(MISSING_VALUE);
        }
    }

    private void CalculateFG(InsilicoMolecule Mol){
        this.setNRCHO(0);this.setNArNO2(0); this.setNRNO2(0);
        try {
            DescriptorBlock FG = new FG();
            FG.Calculate(Mol);
            this.setNRCHO(FG.GetByName("nRCHO").getValue());
            this.setNArNO2(FG.GetByName("nArNO2").getValue());
            this.setNRNO2(FG.GetByName("nRNO2").getValue());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setNRCHO(MISSING_VALUE);this.setNArNO2(MISSING_VALUE); this.setNRNO2(MISSING_VALUE);
        }
    }

    private void CalculateMLogP(InsilicoMolecule Mol){
        this.setMLOGP(0);
        try {
            DescriptorMLogP descriptorMLogP = new DescriptorMLogP(Mol);
            this.setMLOGP(descriptorMLogP.getMLogP());
        } catch (Exception ex){
            log.warn(ex.getMessage());
            this.setMLOGP(MISSING_VALUE);
        }
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

    private boolean isIn(ArrayList<String> list, String s) {
        for (String ss : list)
            if (ss.equalsIgnoreCase(s))
                return true;
        return false;
    }

    /**
     * Sets CATS 2D atom types for each atom, as a list of string containing
     * all matching types for each atom.
     **/
    private ArrayList<String>[] setCatsAtomType(IAtomContainer m) {

        int nSK = m.getAtomCount();
        ArrayList<String>[] AtomTypes = new ArrayList[nSK];

        for (int i=0; i<nSK; i++) {

            AtomTypes[i] = new ArrayList<>();
            IAtom CurAt =  m.getAtom(i);

            boolean tN=false, tP=false, tA=false, tD=false, tL=false;

            // Definition of CATS types
            //
            // A: O, N without H
            // N: [+], NH2
            // P: [-], COOH, POOH, SOOH

            // Hydrogens
            int H = 0;
            try {
                H = CurAt.getImplicitHydrogenCount();
            } catch (Exception e) { }

            // [+]
            if (CurAt.getFormalCharge() > 0) {

                boolean NpOm = false;
                if (ConnAugMatrixCATS[i][i] == 7) {
                    for (int j=0; j<nSK; j++) {
                        if (j==i) continue;
                        if (ConnAugMatrixCATS[i][j]==1) {
                            if (ConnAugMatrixCATS[j][j] == 8) {
                                IAtom Oxy = m.getAtom(j);
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
                tA = true;

                if ( (CurAt.getFormalCharge() == 0) && (H == 1))
                    tD = true;

            }

            // N (NH2 and N without H)
            if (CurAt.getSymbol().equalsIgnoreCase("N")) {

                int nSglBnd = 0, nOtherBnd = 0;
                for (int j=0; j<nSK; j++) {
                    if (j==i) continue;
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1)
                            nSglBnd++;
                        else
                            nOtherBnd++;
                    }
                }

                if ( (CurAt.getFormalCharge() == 0) &&
                        (H == 2) &&
                        (nSglBnd == 1) &&
                        (nOtherBnd == 0) )
                    tP = true;

                if (H == 0)
                    tA = true;

                if  ( (CurAt.getFormalCharge() == 0) &&( (H == 1) || (H ==2) ) )
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if (ConnAugMatrixCATS[i][j] == 1) {
                            nSglBnd++;
                            if (ConnAugMatrixCATS[j][j] == 8) {
                                int Obonds = 0;
                                for (int k=0; k<nSK; k++) {
                                    if (k == j) continue;
                                    if (ConnAugMatrixCATS[k][j]>0) Obonds++;
                                }
                                if (Obonds == 1) nSglOH++;
                            }
                        } else {
                            if ( (ConnAugMatrixCATS[i][j] == 2) && (ConnAugMatrixCATS[j][j] == 8) )
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] > 1.5) ) {
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
                    if (ConnAugMatrixCATS[i][j]>0) {
                        if ( (ConnAugMatrixCATS[j][j] != 6) || (ConnAugMatrixCATS[i][j] != 1) ) {
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


            // Sets final types
            if (tN) AtomTypes[i].add(TYPE_N);
            if (tL) AtomTypes[i].add(TYPE_L);

        }

        return AtomTypes;
    }




}
