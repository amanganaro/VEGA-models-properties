package insilico.mutagenicity_caesar.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.mutagenicity_caesar.descriptors.weights.ACF;
import insilico.mutagenicity_caesar.descriptors.weights.EState;
import insilico.mutagenicity_caesar.descriptors.weights.GhoseCrippenWeights;
import insilico.mutagenicity_caesar.descriptors.weights.VertexDegree;
import lombok.Data;
import lombok.extern.log4j.Log4j;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

@Log4j
@Data
public class EmbeddedDescriptors {

    private double MISSING_VALUE = -999;

    private double Gmin = MISSING_VALUE;
    private double IDWBAR = MISSING_VALUE;
    private double ALogP = MISSING_VALUE;

    public EmbeddedDescriptors(InsilicoMolecule Mol){
        CalculateAllDescriptors(Mol);
    }

    private void CalculateAllDescriptors(InsilicoMolecule Mol) {
        CalculateEStates(Mol);
        CalculateInformationContent(Mol);
        CalculateALogP(Mol);
    }

    private void CalculateALogP(InsilicoMolecule Mol){
        this.setALogP(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setALogP(MISSING_VALUE);
            return;
        }

        // Retrieves or calculates ACF
        DescriptorBlock CurACF  = new ACF();
        CurACF.Calculate(Mol);


        double LogP = 0;
        double[] Frags = CurACF.GetAllValues();

        // Check if some fragments are missing values
        for (double d : Frags)
            if (d == Descriptor.MISSING_VALUE)
                return;

        for (int i=0; i<Frags.length; i++) {
            LogP += Frags[i] * GhoseCrippenWeights.GetHydrophobiticty(i);
        }

        this.setALogP(LogP);


    }


    private void CalculateInformationContent(InsilicoMolecule Mol){
        this.setIDWBAR(0);

        IAtomContainer m;
        int MaxPath = 5;

        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setIDWBAR(MISSING_VALUE);
            return;
        }

        // Gets matrices
        double[][] ConnMat;
        int[][] TopoDistMat;
        try {
            ConnMat = Mol.GetMatrixConnectionAugmented();
            TopoDistMat = Mol.GetMatrixTopologicalDistance();
        } catch (GenericFailureException e) {
            log.warn(e.getMessage());
            this.setIDWBAR(MISSING_VALUE);
            return;
        }

        int nSK = m.getAtomCount();

        int[] VertexDeg = VertexDegree.getWeights(m, true);


        // Information content

        int[] TopoDistFreq = new int[nSK];  // frequencies of topological distances
        double TopoDistFreqSum = 0;
        for (int i=0; i<nSK; i++)
            TopoDistFreq[i] = 0;
        for (int i=0; i<nSK; i++)
            for (int j=i+1; j<nSK; j++) {
                TopoDistFreq[TopoDistMat[i][j]]++;
                TopoDistFreqSum += TopoDistMat[i][j];
            }


        double partial=0;
        for (int j=1; j<nSK; j++)
            partial += (TopoDistFreq[j]) * j * (Math.log(j)/Math.log(2));
        double idw = TopoDistFreqSum * (Math.log(TopoDistFreqSum)/Math.log(2)) - partial;
        double idwbar = idw / TopoDistFreqSum;

        this.setIDWBAR(idwbar);



    }

    private void CalculateEStates(InsilicoMolecule Mol){
        this.setGmin(0);

        IAtomContainer m;
        try {
            m = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            this.setGmin(MISSING_VALUE);
            return;
        }
        int nSK = m.getAtomCount();

        // Get EStates weights
        EState es;
        try {
            es = new EState(Mol.GetStructure());
        } catch (Exception e) {
            this.setGmin(MISSING_VALUE);
            return;
        }

        double Gmin= Descriptor.MISSING_VALUE;

        double Ss = 0;

        for (int at=0; at<m.getAtomCount(); at++) {

            IAtom curAt = m.getAtom(at);

            // Sum of e-states
            Ss += es.getEState()[at];

            // Maximum and minimum Estate/HEstate
            Gmin = (Gmin== MISSING_VALUE) ? es.getEState()[at] : (Math.min(es.getEState()[at], Gmin));
            this.setGmin(Gmin);
        }
    }

}
