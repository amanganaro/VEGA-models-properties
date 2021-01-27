package insilico.daphnia_noec.descriptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * Constitutional descriptors block.<p>
 * NOTE: by default, the molecular weight (MW and AMW) is calculated on
 * carbon-scaled values (carbon mass equal to 1, all other values are scaled
 * on carbon).
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class Constitutional extends DescriptorBlock {

    private static final long serialVersionUID = 1L;
    private static final String BlockName = "Constitutional Descriptors";



    /**
     * Constructor. Sets by default MW calculation with scaled values.
     */
    public Constitutional() {
        super();
        this.Name = Constitutional.BlockName;
    }



    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        this.Add("MW", "");
        this.Add("AMW", "");
        this.Add("Sv", "");
        this.Add("Mv", "");
        this.Add("Sp", "");
        this.Add("Mp", "");
        this.Add("Se", "");
        this.Add("Me", "");

        this.Add("nAt", "");
        this.Add("nSk", "");

        this.Add("nBt", ""); // Total no. of bonds
        this.Add("nBo", ""); // no. of non-H bonds
        this.Add("nBm", ""); // no. of multiple (order>1) bonds
        this.Add("nDblBo", "");
        this.Add("nTrpBo", "");
        this.Add("nArBo", "");
        this.Add("SCBO", "");

        this.Add("nH", "");
        this.Add("nC", "");
        this.Add("nN", "");
        this.Add("nO", "");
        this.Add("nP", "");
        this.Add("nS", "");
        this.Add("nF", "");
        this.Add("nCl", "");
        this.Add("nBr", "");
        this.Add("nI", "");
        this.Add("nB", "");

        this.Add("HPerc", "");
        this.Add("CPerc", "");
        this.Add("NPerc", "");
        this.Add("OPerc", "");
        this.Add("XPerc", "");

        this.Add("nHet", "");
        this.Add("nX", "");

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

        IAtomContainer curMol;
        try {
            curMol = mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            SetAllValues(Descriptor.MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];

            int nTotH=0;
            int nC=0, nN=0, nO=0, nP=0, nS=0;
            int nI=0, nF=0, nCl=0, nBr=0, nB=0;
            int nHet=0;
            double mw=0, amw=0, sv=0, mv=0, sp=0, mp=0, se=0, me=0;


            //// Counts on atoms

            for (int i=0; i<nSK; i++) {

                IAtom CurAt = curMol.getAtom(i);

                // Hydrogens
                H[i] = 0;
                try {
                    H[i] = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                nTotH += H[i];


                if (CurAt.getSymbol().equalsIgnoreCase("C"))
                    nC++;
                else
                    nHet++;

                if (CurAt.getSymbol().equalsIgnoreCase("N"))
                    nN++;
                if (CurAt.getSymbol().equalsIgnoreCase("O"))
                    nO++;
                if (CurAt.getSymbol().equalsIgnoreCase("P"))
                    nP++;
                if (CurAt.getSymbol().equalsIgnoreCase("S"))
                    nS++;
                if (CurAt.getSymbol().equalsIgnoreCase("F"))
                    nF++;
                if (CurAt.getSymbol().equalsIgnoreCase("Cl"))
                    nCl++;
                if (CurAt.getSymbol().equalsIgnoreCase("Br"))
                    nBr++;
                if (CurAt.getSymbol().equalsIgnoreCase("I"))
                    nI++;
                if (CurAt.getSymbol().equalsIgnoreCase("B"))
                    nB++;

            }

            this.SetByName("nAt", nSK + nTotH);
            this.SetByName("nSk", nSK);

            this.SetByName("nH", nTotH);
            this.SetByName("nC", nC);
            this.SetByName("nN", nN);
            this.SetByName("nO", nO);
            this.SetByName("nP", nP);
            this.SetByName("nS", nS);
            this.SetByName("nF", nF);
            this.SetByName("nCl", nCl);
            this.SetByName("nBr", nBr);
            this.SetByName("nI", nI);
            this.SetByName("nB", nB);

            this.SetByName("HPerc", (nTotH/(double)(nSK + nTotH))*100);
            this.SetByName("CPerc", (nC/(double)(nSK + nTotH))*100);
            this.SetByName("NPerc", (nN/(double)(nSK + nTotH))*100);
            this.SetByName("OPerc", (nO/(double)(nSK + nTotH))*100);
            this.SetByName("XPerc", ((nI + nF + nCl + nBr)/(double)(nSK + nTotH))*100);

            this.SetByName("nHet", nHet);
            this.SetByName("nX", nI + nF + nCl + nBr);


            //// Counts on bonds

            int nArBonds=0, nDblBonds=0, nTrpBonds=0, nMulBonds=0;
            double scbo=0;

            for (int i=0; i<nBO; i++) {

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
                        if (CurBo.getOrder() == IBond.Order.TRIPLE) {
                            nTrpBonds++;
                            scbo += 3;
                        }
                    }
                }

            }

            this.SetByName("nBt", nBO + nTotH);
            this.SetByName("nBo", nBO);
            this.SetByName("nBm", nMulBonds);
            this.SetByName("nDblBo", nDblBonds);
            this.SetByName("nTrpBo", nTrpBonds);
            this.SetByName("nArBo", nArBonds);
            this.SetByName("SCBO", scbo);


            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");
            double[] wVdW = VanDerWaals.getWeights(curMol);
            double HVdW = VanDerWaals.GetVdWVolume("H");
            double[] wPol = Polarizability.getWeights(curMol);
            double HPol = Polarizability.GetPolarizability("H");
            double[] wEl = Electronegativity.getWeights(curMol);
            double HEl = Electronegativity.GetElectronegativity("H");

            // for non-scaled molecular weight
//            double[] wMassNS = new double[wMass.length];
//            double CarbonWeight = 12.011;
//            double HMassNS = HMass * CarbonWeight;
//            for (int i=0; i<nSK; i++)
//                wMassNS[i] = wMass[i] * CarbonWeight;

            for (int i=0; i<nSK; i++) {
                if (wMass[i] == -999)
                    mw = -999;
                if (wVdW[i] == -999)
                    sv = -999;
                if (wPol[i] == -999)
                    sp = -999;
                if (wEl[i] == -999)
                    se = -999;
            }

            for (int i=0; i<nSK; i++) {
                if (mw != -999) {
                    mw += wMass[i];
                    if (H[i]>0) {
                        mw += HMass * H[i];
                    }
                }
                if (sv != -999) {
                    sv += wVdW[i];
                    if (H[i]>0)
                        sv += HVdW * H[i];
                }
                if (sp != -999) {
                    sp += wPol[i];
                    if (H[i]>0)
                        sp += HPol * H[i];
                }
                if (se != -999) {
                    se += wEl[i];
                    if (H[i]>0)
                        se += HEl * H[i];
                }
            }

            if (mw != -999)
                amw = mw/(nSK + nTotH);
            if (sv != -999)
                mv = sv/(nSK + nTotH);
            if (sp != -999)
                mp = sp/(nSK + nTotH);
            if (se != -999)
                me = se/(nSK + nTotH);
            this.SetByName("MW", mw);
            this.SetByName("AMW", amw);
            this.SetByName("Sv", sv);
            this.SetByName("Mv", mv);
            this.SetByName("Sp", sp);
            this.SetByName("Mp", mp);
            this.SetByName("Se", se);
            this.SetByName("Me", me);

        } catch (Throwable e) {
            this.SetAllValues(Descriptor.MISSING_VALUE);
        }

    }


    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    public DescriptorBlock CreateClone()
            throws CloneNotSupportedException {
        Constitutional block = new Constitutional();
        block.CloneDetailsFrom(block);
        return block;
    }
}
