package insilico.mutagenicity_amines.rules;

import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;


public class AromaticAminesFunctionalGroups {
    private static final Logger log = LogManager.getLogger(AromaticAminesFunctionalGroups.class);

    private int aggCount;
    private int deactivatingCount;
    private int activatingCount;
    private int conjugatedCount;
    private int narCount;

    private SmartsPattern[] spAgg;
    private SmartsPattern[] spDeactivating;
    private SmartsPattern[] spActivating;
    private SmartsPattern[] spConjugated;
    private SmartsPattern[] spNar;


    private final static String[] smartsAgg = {
            "a[NH2]",
            "a[NH][CH3]",
            "a[N]([CH3])[CH3]",
            "a[NH][CH2][CH3]",
            "a[N]([CH2][CH3])[CH2][CH3]",
            "aN[O,N]",
            "a[NH]C(=O)([CH3])"
    };

    private final static String[] smartsDeactivating= {
            "aC([F,Cl,Br])([F,Cl,Br])[F,Cl,Br]",
            "aC#N",
            "aS(=O)(=O)*"
    };

    private final static String[] smartsActivating = {
            "a[OH]",
            "a[CH3]",
            "a[CH2][CH3]",
            "aO[CH3]",
            "aO[CH2][CH3]"
    };

    private final static String[] smartsConjugated = {
            "aa",
            "aC=Ca"
    };

    private final static String[] smartsNar = {
            "[n]"
    };

    public AromaticAminesFunctionalGroups() {
        aggCount = 0;
        deactivatingCount = 0;
        activatingCount = 0;
        conjugatedCount = 0;
        narCount = 0;

        try {

            spAgg = new SmartsPattern[smartsAgg.length];
            spDeactivating = new SmartsPattern[smartsDeactivating.length];
            spActivating = new SmartsPattern[smartsActivating.length];
            spConjugated = new SmartsPattern[smartsConjugated.length];
            spNar = new SmartsPattern[smartsNar.length];


            int idx = 0;
            for (String smarts : smartsAgg) {
                spAgg[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (spAgg[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : smartsDeactivating) {
                spDeactivating[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (spDeactivating[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : smartsActivating) {
                spActivating[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (spActivating[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : smartsConjugated) {
                spConjugated[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (spConjugated[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : smartsNar) {
                spNar[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (spNar[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

        }  catch (Exception ex){
            log.info(ex.getMessage());
        }

    }

    public void calculateMatches(InsilicoMolecule mol){
        try {
            for(SmartsPattern curPattern: spAgg) {
                if(curPattern.matches(mol.GetStructure()))
                    aggCount++;
            }
            for(SmartsPattern curPattern: spDeactivating) {
                if(curPattern.matches(mol.GetStructure()))
                    deactivatingCount++;
            }
            for(SmartsPattern curPattern: spActivating) {
                if(curPattern.matches(mol.GetStructure()))
                    activatingCount++;
            }
            for(SmartsPattern curPattern: spConjugated) {
                if(curPattern.matches(mol.GetStructure()))
                    conjugatedCount++;
            }
            for(SmartsPattern curPattern: spNar) {
                if(curPattern.matches(mol.GetStructure()))
                    narCount++;
            }
        } catch (InvalidMoleculeException e) {
            log.info(e.getMessage());
        }
    }








    public int getAggCount() {
        return aggCount;
    }

    public int getDeactivatingCount() {
        return deactivatingCount;
    }

    public int getActivatingCount() {
        return activatingCount;
    }

    public int getConjugatedCount() {
        return conjugatedCount;
    }

    public int getNarCount() {
        return narCount;
    }
}
