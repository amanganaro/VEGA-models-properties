package insilico.mutagenicity_amines.rules;

import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.blocks.Constitutional;
import insilico.core.descriptor.blocks.Rings;
import insilico.core.exception.DescriptorNotFoundException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.ringsearch.RingSearch;

import java.util.List;
import java.util.Map;


public class AromaticAminesSubclassClassifier {
    private static final Logger log = LogManager.getLogger(AromaticAminesSubclassClassifier.class);

    public static final String MONOCYCLES = "Monocycles";
    public static final String BICYCLES = "Bicycles, non-fused";
    public static final String THREE_MORE_CYCLES = "Three or more cycles, non-fused";
    public static final String POLYCYCLES_FIVE = "Polycycles, 5 + 6 terms";
    public static final String POLYCYCLES_SIX = "Polycycles, 6 + 6 terms";
    public static final String POLYCYCLES_THREE_FUSED= "Polycycles, more than three fused rings";
    public static final String NOT_CLASSIFICATED = "Not Classificated";

    private String subclassClassification;


    public AromaticAminesSubclassClassifier(InsilicoMolecule mol) {

        try {
            AromaticAtomsCountDescriptor aromaticAtomsCountDescriptor = new AromaticAtomsCountDescriptor();
            int aromaticAtomCount = Integer.parseInt(aromaticAtomsCountDescriptor.calculate(mol.GetStructure()).getValue().toString());

            // The ring bridge count (Rbrid) is the difference between the total ring size
            // and the sum of all the bonds of the ring systems [A.H.Lipkus, J. Chem. Inf. Comput. Sci. 2001, 41, 430-438].

//            IAtomContainer attt = new RingSearch(mol.GetStructure()).ringFragments();
            int ringBridgeCount = Math.abs(new RingSearch(mol.GetStructure()).ringFragments().getAtomCount() - new RingSearch(mol.GetStructure()).ringFragments().getBondCount());

            int numberOfRingSystems = mol.GetSSSR().getAtomContainerCount() - ringBridgeCount;

            if((aromaticAtomCount == 5 || aromaticAtomCount == 6) && numberOfRingSystems == 1 && ringBridgeCount == 0)
                setSubclassClassification(MONOCYCLES);
            else if((aromaticAtomCount == 10 || aromaticAtomCount == 11 || aromaticAtomCount == 12) && numberOfRingSystems == 2 && ringBridgeCount == 0)
                setSubclassClassification(BICYCLES);
            else if(numberOfRingSystems >=3 && ringBridgeCount == 0)
                setSubclassClassification(THREE_MORE_CYCLES);
            else if(aromaticAtomCount == 9 && numberOfRingSystems == 1 && ringBridgeCount == 1)
                setSubclassClassification(POLYCYCLES_FIVE);
            else if(aromaticAtomCount == 10 && numberOfRingSystems == 1 && ringBridgeCount == 1)
                setSubclassClassification(POLYCYCLES_SIX);
            else if(numberOfRingSystems == 1 && ringBridgeCount >= 2)
                setSubclassClassification(POLYCYCLES_THREE_FUSED);

            else setSubclassClassification(NOT_CLASSIFICATED);

        } catch (InvalidMoleculeException ex){
            log.warn(ex.getClass() + ": " + ex.getMessage());
            setSubclassClassification(NOT_CLASSIFICATED);
        }
    }

    public String getClassification() {
        return subclassClassification;
    }

    public void setSubclassClassification(String subclassClassification){
        this.subclassClassification = subclassClassification;
    }
}
