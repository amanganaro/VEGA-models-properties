package insilico.mutagenicity_amines.rules;

import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.log4j.Log4j;

@Log4j
public class AromaticAminesModel {

    private final AromaticAminesDescriptors aromaticAminesDescriptors;
    private final AromaticAminesSubclassClassifier aromaticAminesSubclassClassifier;
    private final AromaticAminesFunctionalGroups aromaticAminesFunctionalGroups;

    public static final short NON_MUTAGENIC = 0;
    public static final short MUTAGENIC = 1;
    public static final short NA = -1;



    public AromaticAminesModel(InsilicoMolecule mol){

        aromaticAminesDescriptors = new AromaticAminesDescriptors(mol);
        aromaticAminesSubclassClassifier = new AromaticAminesSubclassClassifier(mol);
        aromaticAminesFunctionalGroups = new AromaticAminesFunctionalGroups();

        aromaticAminesFunctionalGroups.calculateMatches(mol);
    }

    public short calculatePrediction() throws Exception {

        // first condition block
        if(aromaticAminesDescriptors.getDescriptorBasedClassification())
            return NON_MUTAGENIC;

        else {
            // classification into subclasses
            switch (aromaticAminesSubclassClassifier.getClassification()){

                case (AromaticAminesSubclassClassifier.MONOCYCLES):
                    return calculateMonocyclesPrediction();
                case(AromaticAminesSubclassClassifier.BICYCLES):
                    return calculateBicyclesPrediction();
                case(AromaticAminesSubclassClassifier.THREE_MORE_CYCLES):
                    if(aromaticAminesDescriptors.getMw() > 320)
                        return NON_MUTAGENIC;
                    else return MUTAGENIC;
                case (AromaticAminesSubclassClassifier.POLYCYCLES_FIVE):
                case (AromaticAminesSubclassClassifier.POLYCYCLES_THREE_FUSED):
                    return MUTAGENIC;
                case (AromaticAminesSubclassClassifier.POLYCYCLES_SIX):
                    if (aromaticAminesFunctionalGroups.getDeactivatingCount() > 0)
                        return NON_MUTAGENIC;
                    else return MUTAGENIC;
                default:
                    return NA;
            }
        }

    }

    private short calculateMonocyclesPrediction() {


        if(aromaticAminesFunctionalGroups.getAggCount() > 0){
            if(aromaticAminesFunctionalGroups.getNarCount() >= 1)
                return NON_MUTAGENIC;
            else return MUTAGENIC;
        } else {
            if(aromaticAminesFunctionalGroups.getNarCount() >= 1 || aromaticAminesFunctionalGroups.getDeactivatingCount() > 0)
                return NON_MUTAGENIC;
            else {
                if(aromaticAminesFunctionalGroups.getActivatingCount() > 0)
                    return MUTAGENIC;
                else return NON_MUTAGENIC;
            }
        }
    }

    private short calculateBicyclesPrediction() {

        // if have conjugated structures
        if(aromaticAminesFunctionalGroups.getConjugatedCount() > 0){
            if((aromaticAminesFunctionalGroups.getDeactivatingCount() > aromaticAminesFunctionalGroups.getAggCount()) || ((aromaticAminesFunctionalGroups.getNarCount() + 1) > aromaticAminesFunctionalGroups.getAggCount()))
                return NON_MUTAGENIC;
            else return MUTAGENIC;
        } else {
            if(aromaticAminesFunctionalGroups.getAggCount() > 0)
                return MUTAGENIC;
            else return NON_MUTAGENIC;
        }
    }


}
