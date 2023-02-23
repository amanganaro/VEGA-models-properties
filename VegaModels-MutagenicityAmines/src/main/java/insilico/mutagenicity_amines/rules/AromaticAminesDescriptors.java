package insilico.mutagenicity_amines.rules;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.tools.AtomicNumber;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.beam.Element;


public class AromaticAminesDescriptors {
    private static final Logger log = LogManager.getLogger(AromaticAminesDescriptors.class);

    private double MISSING_VALUE = -999;

    private double Mw = MISSING_VALUE;
    private double nRotableBonds = MISSING_VALUE;
    private double largestChain = MISSING_VALUE;
    private double sp3Ratio = MISSING_VALUE;

    public AromaticAminesDescriptors(InsilicoMolecule mol){

        try {
            DescriptorBlock constitutional = new Constitutional();
            constitutional.Calculate(mol);
            Mw = constitutional.GetByName("MW_da").getValue();

        } catch (DescriptorNotFoundException e) {
            log.warn(e.getMessage());
        }

        try {
            RotatableBondsCountDescriptor rotatableBondsCountDescriptor = new RotatableBondsCountDescriptor();
            nRotableBonds = Double.parseDouble(rotatableBondsCountDescriptor.calculate(mol.GetStructure()).getValue().toString());
        } catch (Exception e){
            log.warn(e.getClass() + " - " + e.getMessage());
        }
        try {
            LargestChainDescriptor largestChainDescriptor = new LargestChainDescriptor();
            largestChain = Double.parseDouble(largestChainDescriptor.calculate(mol.GetStructure()).getValue().toString());
        } catch (Exception e){
            log.warn(e.getClass() + " - " + e.getMessage());
        }

        try {
            sp3Ratio = calculateSp3Ratio(mol.GetStructure());
        } catch (InvalidMoleculeException e) {
            log.warn(e.getClass() + " - " + e.getMessage());
        }

    }

    public boolean getDescriptorBasedClassification() throws DescriptorNotFoundException {
        if(Mw != MISSING_VALUE && nRotableBonds != MISSING_VALUE && largestChain != MISSING_VALUE && sp3Ratio != MISSING_VALUE)
            return (Mw > 430 && nRotableBonds > 5 && largestChain > 10 && sp3Ratio > 0.17);
        else throw new DescriptorNotFoundException("One or more descriptors are missing");
    }

    private double calculateSp3Ratio(IAtomContainer container){
        try {
            IAtomContainer clone = container.clone();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clone);
            int nC = 0;
            int nsp3 = 0;
            for (IAtom atom : clone.atoms()) {
                if (atom.getAtomicNumber() != Element.Carbon.atomicNumber())
                    continue;
                else {
                    nC++;
                    if (atom.getHybridization() == IAtomType.Hybridization.SP3) nsp3++;
                }
            }

            if (nC+nsp3 == 0) {
                return 0;
            }
            return nsp3 / (double) (nC);
        } catch (CloneNotSupportedException | CDKException e) {
            log.warn(e.getClass() + ": " + e.getMessage());
            return MISSING_VALUE;
        }
    }

    public double getMw() {
        return Mw;
    }

    public void setMw(double mw) {
        Mw = mw;
    }

    public double getnRotableBonds() {
        return nRotableBonds;
    }

    public void setnRotableBonds(double nRotableBonds) {
        this.nRotableBonds = nRotableBonds;
    }

    public double getLargestChain() {
        return largestChain;
    }

    public void setLargestChain(double largestChain) {
        this.largestChain = largestChain;
    }

    public double getSp3Ratio() {
        return sp3Ratio;
    }

    public void setSp3Ratio(double sp3Ratio) {
        this.sp3Ratio = sp3Ratio;
    }
}
