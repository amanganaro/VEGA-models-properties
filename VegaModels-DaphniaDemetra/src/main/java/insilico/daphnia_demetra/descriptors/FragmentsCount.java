package insilico.daphnia_demetra.descriptors;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.ArrayList;
import java.util.List;

public class FragmentsCount {

    private Pattern q1;

    public int O_057;  //  Phenol/Enol/Carboxyl OH and P-OH groups
    public int O_060;  // da definire meglio!
    public int S_107;  // R2S / RS-SR
    public int nHAcc;  // number of N, O, F
    public int nArNR2; // number of nArNR2 groups

    public FragmentsCount() {
        O_057 = -999;
        O_060 = -999;
        S_107 = -999;
        nHAcc = -999;
        nArNR2 = -999;
        try {
            q1 = SmartsPattern.create("[a][N;D3]([$([C,c]);!$([C,c]=O)])[$([C,c]);!$([C,c]=O)]", DefaultChemObjectBuilder.getInstance());
            ((SmartsPattern)q1).setPrepare(false);
        } catch (Exception e){
            q1 = null;
        }
    }

    private static int CountHydrogens(IAtomContainer mol, int atom) {

        int H=0;
        if (mol.getAtom(atom).getImplicitHydrogenCount()!=null)
            H = mol.getAtom(atom).getImplicitHydrogenCount();
        return H;
    }

    public void Calculate(InsilicoMolecule Mol) throws GenericFailureException {

        IAtomContainer mol;
        try {
            mol = Mol.GetStructure();
        } catch (InvalidMoleculeException ex){
            throw new GenericFailureException(ex);
        }

        O_057 = 0;
        O_060 = 0;
        S_107 = 0;
        nHAcc = 0;
        nArNR2 = 0;

        int n = mol.getAtomCount();
        double[][] ConnMatrix = Mol.GetMatrixConnectionAugmented();

        RingSet MolRings;
        try {
            MolRings = Mol.GetSSSR();
        } catch (InvalidMoleculeException e) {
            throw new GenericFailureException(e.getMessage());
        }

        try {

            if (q1 == null)
                throw new Exception("");
            boolean status = q1.matches(Mol.GetStructure());
            int nmatch = 0;
            List<Mappings> mappings = new ArrayList<>();
            if (status) {
                mappings.add(q1.matchAll(Mol.GetStructure()));
                nmatch = q1.matchAll(Mol.GetStructure()).countUnique();
            }
            nArNR2 = nmatch;

        } catch (Exception e) {
            throw new GenericFailureException(e.getMessage());
        }

        ArrayList VisitedS = new ArrayList();

        for(int i = 0; i < n ; i++) {

            if(mol.getAtom(i).getSymbol().equalsIgnoreCase("O")) {
                nHAcc++;

                if (CountHydrogens(mol, i) == 1) {
                    //// -OH Group
                    int ConnAtomIdx = -1;
                    for (int j = 0; j < n; j++) {
                        if (j == i) continue;
                        if (ConnMatrix[i][j] > 0) {
                            ConnAtomIdx = j;
                            break;
                        }
                    }
                    if (mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("C")) {

                        int Carom = 0, Csingle = 0, Cdouble = 0, Odouble = 0, Other = 0;
                        for (int j = 0; j < n; j++) {
                            if ((j == ConnAtomIdx) || (j == i)) continue;
                            if (ConnMatrix[ConnAtomIdx][j] > 0) {
                                if (mol.getAtom(j).getSymbol().equalsIgnoreCase("C")) {
                                    if (ConnMatrix[ConnAtomIdx][j] == 1) Csingle++;
                                    if (ConnMatrix[ConnAtomIdx][j] == 2) Cdouble++;
                                    if (ConnMatrix[ConnAtomIdx][j] == 1.5) Carom++;
                                } else {
                                    if ((mol.getAtom(j).getSymbol().equalsIgnoreCase("O")) &&
                                            (ConnMatrix[ConnAtomIdx][j] == 2))
                                        Odouble++;
                                    else
                                        Other++;
                                }
                            }
                        }

                        // phenol / enol / carboxyl -OH
                        if (((Carom == 2) && (Other == 0)) || ((Csingle == 1) && (Cdouble == 1) && (Other == 0)) || ((Odouble == 1) && ((Other == 1) || (Csingle == 1))))
                            O_057++;

                    } else if (mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("P")) {
                        O_057++;
                    }
                } else if (CountHydrogens(mol, i) == 0) {

                    ////
                    //// Possible -O- Group
                    ////

                    // O in ring
                    // if ..O.. increase counter otherwise ignore
                    if (MolRings.contains(mol.getAtom(i)))
                        if (mol.getAtom(i).getFlag(CDKConstants.ISAROMATIC)) {
                            O_060++;
                            continue;
                        }

                    int ConnAtomIdx;
                    int ConnCount = 0, Carom = 0, CdoubleX = 0, P_S = 0;
                    for (int j = 0; j < n; j++) {
                        if (j == i) continue;
                        if (ConnMatrix[i][j] > 0) {
                            ConnAtomIdx = j;
                            ConnCount++;

                            // Checks what the connected atom is
                            if (mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("C")) {

                                // Aromatic Carbon
                                if (mol.getAtom(ConnAtomIdx).getFlag(CDKConstants.ISAROMATIC)) {
                                    Carom++;
                                } else {
                                    for (int k = 0; k < n; k++) {
                                        if ((k == ConnAtomIdx) || (k == i)) continue;
                                        if (ConnMatrix[ConnAtomIdx][k] > 0)
                                            if (ConnMatrix[ConnAtomIdx][k] == 2) {

                                                // -O-C=X fragment
                                                // Checks if C=X is in the same ring
                                                if (mol.getAtom(k).getSymbol().equalsIgnoreCase("C")) {
                                                    RingSet rs = (RingSet) MolRings.getRings(mol.getAtom(ConnAtomIdx));
                                                    if (!(rs.contains(mol.getAtom(k))))
                                                        CdoubleX++;
                                                } else
                                                    CdoubleX++;
                                            }
                                    }
                                }
                            } else {

                                if ((mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("P")) ||
                                        (mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("S")))
                                    P_S++;
                            }
                        }
                    }

                    // Checks for O_60 fragment:
                    if (ConnCount == 2) {

                        // -O-P
                        if (P_S == 1) {
                            O_060++;
                            continue;
                        }

                        // Ar-O-Ar
                        if (Carom == 2) {
                            O_060++;
                            continue;
                        }

                        // -O-C=X
                        if ((CdoubleX == 1)) {
                            O_060++;
                        }
                    }
                }
            }
            if (mol.getAtom(i).getSymbol().equalsIgnoreCase("S")) {

                ////
                //// -S- Group
                ////

                int ConnCount=0, R=0, S_SR=0;

                VisitedS.add(i);

                for (int ConnAtomIdx=0; ConnAtomIdx<n; ConnAtomIdx++){
                    if (i==ConnAtomIdx) continue;
                    if (ConnMatrix[i][ConnAtomIdx]>0) {
                        ConnCount++;
                        if (mol.getAtom(ConnAtomIdx).getSymbol().equalsIgnoreCase("S")) {
                            boolean AlreadyVisited = false;
                            for (int w=0; w<VisitedS.size(); w++)
                                if ((Integer) VisitedS.get(w) ==ConnAtomIdx) {
                                    AlreadyVisited = true;
                                    break;
                                }
                            if (!(AlreadyVisited)) {
                                VisitedS.add(ConnAtomIdx);
                                for (int k=0; k<n; k++){
                                    if (k==ConnAtomIdx) continue;
                                    if (ConnMatrix[k][ConnAtomIdx]>0) {
                                        S_SR++;
                                        break;
                                    }
                                }
                            }
                        } else {
                            R++;
                        }

                    }
                }

                if (ConnCount==2) {
                    // R2S
                    if (R==2) {
                        S_107++;
                        continue;
                    }

                    // RS-SR
                    if (S_SR==1) {
                        S_107++;
                        continue;
                    }
                }

                continue;
            }


            if (mol.getAtom(i).getSymbol().equalsIgnoreCase("N")) {
                nHAcc++;
                continue;
            }


            if (mol.getAtom(i).getSymbol().equalsIgnoreCase("F")) {
                nHAcc++;
            }
        }
    }

}
