package insilico.logk.descriptors.weights;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.matrix.TopoDistanceMatrix;
import insilico.core.tools.utils.MoleculeUtilities;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtomContainer;

@Slf4j
public class MoleculePaths {

    // Limit for the number of atoms
    private static final int MAX_PATH_LENGTH = 2000;
    private static final int MAX_PATH_LENGTH_FOR_WALKS = 10;

    private IAtomContainer m;

    // Walk indices
    public double[] Walk_Counts;
    public double[] Self_Returning_Walk_Counts;
    public double Total_Walk_Count;

    // Path indices
    public double[] Path_Counts;
    public double[] Multiple_Path_Counts;
    public double[] TotalPC;
    public double Total_Path_Count;
    public double IDpi;
    public double PCR;
    public double PCD;
    public double ID_Randic;
    public double ID_Balaban;
    public double[] Pws;


    // private vars accessible to all methods for iterative DFS
    private double[] Vertex_Distance_Degree;
    private boolean[] Entered;
    private double[] TotPC;
    private double[] TotPCMult;
    private double[] IDRandic;
    private double[] IDBalaban;
    private double[][] AdjConnectionMatrix;
    private int nSK;



    /**
     * Constructor. When the object is created it directly calculates all path indices.
     * The input molecular structure should be provided WITHOUT explicit hydrogens.
     *
     * @param Mol
     * @throws GenericFailureException
     */
    public MoleculePaths(InsilicoMolecule Mol) throws GenericFailureException {

        Calculate(Mol);

    }


    private void Calculate(InsilicoMolecule Mol) throws GenericFailureException {

        //// Get structure, matrices and initial settings

        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            log.warn("Invalid structure, unable to calculate paths");
            throw new GenericFailureException("invalid structure");
        }

        int[][] AdjMat;
        try {
            AdjMat = Mol.GetMatrixAdjacency();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            throw new GenericFailureException(e.getMessage());
        }
        double[][] AdjMatDbl = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjMatDbl[i][j] = AdjMat[i][j];
        double[][] AdjxAdj = new double[AdjMat.length][AdjMat[0].length];
        for (int i=0; i<AdjMat.length; i++)
            for (int j=0; j<AdjMat[0].length; j++)
                AdjxAdj[i][j] = AdjMat[i][j];

        int[][] TopoDistMatrix;
        try {
            TopoDistMatrix = TopoDistanceMatrix.getMatrix(m);
        } catch (Exception e) {
            log.warn(e.getMessage());
            throw new GenericFailureException(e.getMessage());
        }
        try {
            AdjConnectionMatrix = Mol.GetMatrixConnectionAugmented();
        } catch (Exception e) {
            log.warn(e.getMessage());
            throw new GenericFailureException(e.getMessage());
        }



        nSK = m.getAtomCount();
        int nBO = m.getBondCount();

        Walk_Counts = new double[MAX_PATH_LENGTH];
        Self_Returning_Walk_Counts = new double[MAX_PATH_LENGTH];
        Path_Counts = new double[MAX_PATH_LENGTH];
        Multiple_Path_Counts = new double[MAX_PATH_LENGTH];

        double[][] AdjPow = new double[nSK][nSK];

        double[][] AWC = new double[nSK][MAX_PATH_LENGTH];


        //// Walk calculation + AWC (needed in paths)

        Total_Walk_Count = nSK + nBO;
        Walk_Counts[0] = nBO;
        Self_Returning_Walk_Counts[0] = nSK;

        int Pow_k = 1;
        while (Pow_k < MAX_PATH_LENGTH_FOR_WALKS) {

            for (int i=0; i<nSK; i++) {
                for (int j = 0; j < nSK; j++) {
                    double ww = 0;
                    for (int m = 0; m < nSK; m++)
                        ww += AdjMatDbl[i][m] * AdjxAdj[m][j];
                    AdjPow[i][j] = ww;
                }
            }

            int ww = 0;
            for (int i=0; i<nSK; i++) {
                Self_Returning_Walk_Counts[Pow_k] = Self_Returning_Walk_Counts[Pow_k] + AdjPow[i][i];
                int aw = 0;
                for (int j=0; j< nSK; j++)
                    aw += AdjPow[i][j];
                AWC[i][Pow_k] = aw;
                ww += aw;
            }

            Walk_Counts[Pow_k] = ww;
            Self_Returning_Walk_Counts[Pow_k] = Self_Returning_Walk_Counts[Pow_k];

            Total_Walk_Count = Total_Walk_Count + Walk_Counts[Pow_k];

            Pow_k++;

            for (int i=0; i<nSK; i++)
                for (int j = 0; j < nSK; j++)
                    AdjxAdj[i][j] = AdjPow[i][j];

        }


        //// Paths calculation

        int LenChi = 5;

        int lim = Math.min(nSK, MAX_PATH_LENGTH);
        TotPC = new double[lim];
        TotalPC = new double[lim];
        TotPCMult = new double[lim];
        IDRandic = new double[lim];
        IDBalaban = new double[lim];
        for (int z=0; z<lim;z++)
            TotPC[z] = 0;

        Vertex_Distance_Degree = new double[nSK];

        for (int i=0; i<nSK; i++) {
            Vertex_Distance_Degree[i] = 0;
            for (int j = 0; j < nSK; j++)
                Vertex_Distance_Degree[i] += (double) TopoDistMatrix[i][j];
        }

        Entered = new boolean[nSK];
        int jlim = LenChi;
        Pws = new double[jlim];
        for (int i=0; i<jlim; i++)
            Pws[i] = 0;
        if (nSK <= LenChi)
            jlim = nSK - 1;
        double[] TotPCprev = new double[jlim];

        int ii = 0;

        for (int i=0; i<nSK; i++) {

            for (int j = 0; j < nSK; j++)
                Entered[j] = false;
            Entered[i] = true;

            int PathLength = 0;

            double VD1 = 0;
            for (int k=0; k<nSK; k++)
                if (k != i)
                    if (AdjMat[i][k] != 0)
                        VD1++;

            // Cycle on all atoms connected to i-th
            for (int j = 0; j < nSK; j++) {
                if (j == i) continue;
                if (AdjConnectionMatrix[i][j] != 0) {

                    double VD2 = 0;
                    for (int k=0; k<nSK; k++)
                        if (k != j)
                            if (AdjMat[j][k] != 0)
                                VD2++;

                    double Cur_Mult_Bond_Order = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(i), m.getAtom(j)));
                    double Cur_Mult_Ver_Deg = 1.0 / Math.sqrt(VD1 * VD2);
                    double Cur_Balaban_Weight = 1.0 / Math.sqrt(Vertex_Distance_Degree[i] * Vertex_Distance_Degree[j]);

                    NextPathVisit(j, PathLength, Cur_Mult_Bond_Order, Cur_Mult_Ver_Deg, Cur_Balaban_Weight);

                }
            }

            for (int k=0; k<jlim; k++) {
                if (AWC[ii][k] != 0)
                    Pws[k] += (TotPC[k] - TotPCprev[k]) / AWC[ii][k];
                TotPCprev[k] = TotPC[k];
            }

            ii++;
        }

        for (int k=0; k<jlim; k++)
            Pws[k] = Pws[k] / nSK;

        for (int i=0; i<lim; i++) {
            TotPC[i] /= 2.0;
            TotalPC[i] = TotPC[i];
            TotPCMult[i] /= 2.0;
        }


        Total_Path_Count = nSK;
        IDpi = nSK;
        ID_Randic = nSK;
        ID_Balaban = nSK;

        for (int i=0; i<lim; i++) {
            Total_Path_Count = Total_Path_Count + TotPC[i];
            IDpi = IDpi + TotPCMult[i];
            ID_Randic = ID_Randic + (IDRandic[i] / 2.0);
            ID_Balaban = ID_Balaban + (IDBalaban[i] / 2.0);
        }

        PCR = (IDpi / Total_Path_Count);
        PCD = (IDpi - Total_Path_Count);

        for (int i=0; i<lim; i++) {
            Path_Counts[i] = TotPC[i];
            Multiple_Path_Counts[i] = TotPCMult[i];
        }
    }


    private void NextPathVisit(int Atom_Idx, int PathLength, double Mult_Bond_Order,
                               double Mult_Ver_Deg, double Balaban_Weight) {

        Entered[Atom_Idx] = true;
        if (PathLength < TotPC.length) {

            TotPC[PathLength] += 1;
            TotPCMult[PathLength] += Mult_Bond_Order;
            IDRandic[PathLength] += Mult_Ver_Deg;
            IDBalaban[PathLength] += Balaban_Weight;

            for (int Next_Atom_Idx = 0; Next_Atom_Idx < nSK; Next_Atom_Idx++) {
                if (Next_Atom_Idx == Atom_Idx) continue;
                if (Entered[Next_Atom_Idx]) continue;
                if (AdjConnectionMatrix[Atom_Idx][Next_Atom_Idx] != 0) {

                    double BW1 = Vertex_Distance_Degree[Atom_Idx];
                    double BW2 = Vertex_Distance_Degree[Next_Atom_Idx];

                    double BO = MoleculeUtilities.Bond2Double(m.getBond(m.getAtom(Atom_Idx), m.getAtom(Next_Atom_Idx)));

                    double VD = 0;
                    for (int k=0; k<nSK; k++)
                        if (k != Atom_Idx)
                            if (AdjConnectionMatrix[Atom_Idx][k] != 0)
                                VD++;

                    double VD2= 0;
                    for (int k=0; k<nSK; k++)
                        if (k != Next_Atom_Idx)
                            if (AdjConnectionMatrix[Next_Atom_Idx][k] != 0)
                                VD2++;



                    double Cur_Mult_Bond_Order = BO * Mult_Bond_Order;
                    double Cur_Mult_Ver_Deg = Mult_Ver_Deg * 1.0 / Math.sqrt(VD * VD2);
                    double Cur_Balaban_Weight = Balaban_Weight * 1.0 / Math.sqrt(BW1 * BW2);

                    NextPathVisit(Next_Atom_Idx, PathLength + 1, Cur_Mult_Bond_Order, Cur_Mult_Ver_Deg, Cur_Balaban_Weight);

                }
            }

        }
        Entered[Atom_Idx] = false;

    }

}
