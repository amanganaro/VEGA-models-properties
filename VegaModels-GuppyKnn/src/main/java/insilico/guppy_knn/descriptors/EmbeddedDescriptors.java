package insilico.guppy_knn.descriptors;

import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.guppy_knn.descriptors.weights.Electronegativity;
import insilico.guppy_knn.descriptors.weights.Mass;
import insilico.guppy_knn.descriptors.weights.Polarizability;
import insilico.guppy_knn.descriptors.weights.VanDerWaals;
//import insilico.logp_alogp.descriptors.weights.ACF;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Arrays;
import java.util.List;

@Slf4j
@Data
public class EmbeddedDescriptors {

    public final int MISSING_VALUE = -999;

    private double MW = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol) {
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) {
        CalculateMW(Mol);
    }

    private void CalculateMW(InsilicoMolecule Mol){
        this.setMW(0);
        IAtomContainer curMol;
        try {
            curMol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setMW(MISSING_VALUE);
            return;
        }

        try {

            int nSK = curMol.getAtomCount();
            int nBO = curMol.getBondCount();
            int[] H = new int[nSK];


            double mw=0, amw=0, sv=0, mv=0, sp=0, mp=0, se=0, me=0;

            // Weights sums
            double[] wMass = Mass.getWeights(curMol);
            double HMass = Mass.GetMass("H");
            double[] wVdW = VanDerWaals.getWeights(curMol);
            double HVdW = VanDerWaals.GetVdWVolume("H");
            double[] wPol = Polarizability.getWeights(curMol);
            double HPol = Polarizability.GetPolarizability("H");
            double[] wEl = Electronegativity.getWeights(curMol);
            double HEl = Electronegativity.GetElectronegativity("H");



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
            }

            this.setMW(mw);


        } catch (Throwable e) {
            this.setMW(MISSING_VALUE);
        }
    }


}
