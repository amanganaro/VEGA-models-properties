package insilico.glucocorticoid_receptor.alert;

import insilico.core.alerts.AlertBlockFromSMARTS;
import insilico.core.alerts.AlertList;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

@Slf4j
public class SAGlucocorticoidReceptor {

    private final static short INACTIVE = 0;
    private final static short ACTIVE_AGONIST = 1;
    private final static short ACTIVE_ANTAGONIST = 2;
    private final static short ACTIVE_A_ANTA = 3;
    private final static short NOT_PREDICTED = -1;

    private Short Matches;
    private List<String> SAMatched;

    private SmartsPattern[] SA_SmartsPrimary_I;
    private SmartsPattern[] SA_SmartsPrimary_II;
    private SmartsPattern[] SA_SmartsPrimary_III;
    private SmartsPattern[] SA_SmartsSecondary_1_1;
    private SmartsPattern[] SA_SmartsSecondary_2_1;
    private SmartsPattern[] SA_SmartsSecondary_3_1;
    private SmartsPattern[] SA_SmartsTertiary_1_1_1;
    private SmartsPattern[] SA_SmartsTertiary_1_1_2;
    private SmartsPattern[] SA_SmartsTertiary_1_1_3;
    private SmartsPattern[] SA_SmartsTertiary_1_1_4;
    private SmartsPattern[] SA_SmartsTertiary_2_1_1;
    private SmartsPattern[] SA_SmartsTertiary_3_1_1;
    private SmartsPattern[] SA_SmartsTertiary_3_1_1_1;
    private SmartsPattern[] SA_SmartsTertiary_3_2_1;
    private SmartsPattern[] SA_SmartsTertiary_3_2_1_1;
    private SmartsPattern[] SA_SmartsTertiary_3_3_1;

//
//    C-N-C:C
//    N-C:C-[#1]
//    N-C:C-C
//    N-C-C:C
//    C-N-C-C:C
//    N-C-C-C:C
//    N-C=C-[#1]
//    N-C-C=C
//    C-C-N-C-C
//    C-N-C-C-C
//    N-C-C-C-C-C
//    C-C(C)-C(C)-C
//    C-C-C-C-C-C(C)-C

//    private final static String[] SmartsPrimaryI = {
//            "C-N-[C,c]:[C,c]",
//            "N-[C,c]:[CH,cH]",
//            "N-[C,c]:[C,c]-C",
//            "N-C-[C,c]:[C,c]",
//            "C-N-C-[C,c]:[C,c]",
//            "N-C-C-[C,c]:[C,c]",
//    };
    private final static String[] SmartsPrimaryI = {
            "[#6][#7][#6]:[#6]",
            "[#7][#6]:[CH,cH]",
            "[#7][#6]:[#6][#6]",
            "[#7][#6][#6]:[#6]",
            "[#6][#7][#6][#6]:[#6]",
            "[#7][#6][#6][#6]:[#6]",
    };

//    private final static String[] SmartsPrimaryII = {
////            "N-C=[CH,cH]",
//            "N-C=[CH,cH]",
//            "N-C-C=C",
//            "C-C-N-C-C",
//            "C-N-C-C-C",
//            "N-C-C-C-C-C",
//    };

    private final static String[] SmartsPrimaryII = {
//            "N-C=[CH,cH]",
            "[#7][#6]=[CH,cH]",
            "[#7][#6][#6]=[#6]",
            "[#6][#6][#7][#6][#6]",
            "[#6][#7][#6][#6][#6]",
            "[#7][#6][#6][#6][#6][#6]",
    };

//    private final static String[] SmartsPrimaryIII = {
//            "C-C(C)-C(C)-C",
//            "C-C-C-C-C-C(C)-C"
//    };

    private final static String[] SmartsPrimaryIII = {
            "[#6][#6]([#6])[#6]([#6])[#6]",
            "[#6][#6][#6][#6][#6][#6]([#6])[#6]"
    };

    private final static String[] SmartsSecondary_1_1 = {
            "C(C)(C)CC(C)(O)",
            "CCc1cc(ccc1)c1ccccc1",
            "c1cnn(c1)c1ccc(F)cc1",
            "CC(C)(CN)C(c1ccccc1)",
            "CC(C)NS(=O)(=O)",
            "c1cccc(NS(=O)(=O)C)c1",
            "CC1CCNc2ccc(cc12)",
            "O=c1[nH]c(=O)c(c([nH]1))C",
            "CC(CC(O))c1ccccc1",
            "COc1cccc2OCc3cc(N)ccc3c12",
            "C(Nc1ccc(cc1))c1cccc(F)c1",
            "CC1=CC(C)(C)Nc2ccc(cc12)"
    };

    private final static String[] SmartsSecondary_2_1 = {
            "O=C1NC=C(Cc2ccccc2)C(=O)N1",
            // "n1c(=O)[nH]c(c(Cc2ccccc2)c1=O)",  this has been changed with the prev one, structure can not be aromatic
            "CC(C)(C(c1ccccc1))CN"
    };

    private final static String[] SmartsSecondary_3_1 = {
            "C12CCC(CC1CCc1cc(O)ccc21)",
            "C1CCCC2=CC(=O)CCC12",
            "C(C)C(F)",
            "CCC(C)(c1ccc(O)cc1)",
            "CC2C(=C)CCCC2",
            "CC2CCCCC2(C)C",
            "CC12CCC(CC1CCc1cc(O)ccc21)",
            "C(O)C1CCCC2=CC(=O)CCC12C",
            "SCCO"
    };

    private final static String[] SmartsTertiary_1_1_1 = {
            "CC(C)NS(=O)(=O)"
    };

    private final static String[] SmartsTertiary_1_1_2 = {
            "c1cccc(NS(=O)(=O)C)c1",
            "CC1=CNC(=O)NC1=O",
            // "O=c1[nH]c(=O)c(c([nH]1))C",  this has been changed with the prev one, structure can not be aromatic
            "C(Nc1ccc(cc1))c1cccc(F)c1"
    };

    private final static String[] SmartsTertiary_1_1_3 = {
            "CC(C)C(=O)NC(=O)",
            "CC(O)(Cc1cc2cc(N)ncc2[nH]1)C(F)(F)F",
            "Cc1cc2cc(ncc2[nH]1)S(=O)(=O)C",
            "CC1C(O)C(C)(C)Nc2c(C)cc(c(C)c12)",
            "c1ccc(F)cc1C(C)(C)CC(O)C(C)C",
            "C(O)(CNC(=O)c1c(F)cccc1)C(F)(F)F",
            "c1ccc2c(cnn2C)c1",
            "CC(O)(Cc1cc2ccc(C)cc2[nH]1)C(F)(F)F",
            "CN(CC(O)C(F)(F)F)C(=O)",
            "CC(O)(Cc1cc2ccccc2[nH]1)C(F)(F)F",
            "OC(c1ccccc1)C(F)(F)F",
            "CCCC(C)(C)C(=O)N",
            "CC1C(O)C(C)(C)Nc2c(C)cc(cc12)",
            "Cc1ccc2c(cnn2)c1",
            "C(O)(CNC(=O))C(F)(F)F",
            "Cc3sccc3C",
            "CC(C)CCC(C)C(=O)",
            "C1(CC(O)C(F)(F)F)CCCc2ccccc12",
            "CN(c1ccccc1)S(=O)(=O)",
            "Cc1cc2cc(CC(O)(C)C(F)(F)F)[nH]c2cn1",
            "CC1C(O)C(C)(C)Nc2c(F)cc(cc12)"
    };

    private final static String[] SmartsTertiary_1_1_4 = {
            "CC1C(O)C(C)(C)Nc2cc(F)c(c(F)c12)",
            "Cc1ccccc1CN",
            "CN(C)C(=O)c1ccc(cc1)",
            "CC(C)(C(c1ccccc1)c1ccc(O)cc1)",
            "c1ccc2c(ncn2c1)c1ccccc1",
            "CC(C)(C)c3cccc4CCOc34",
            "CC(C)(C(c1ccccc1))C(=O)Nc1nncs1",
            "CC(C)c1cc(ccc1)c1ccccc1",
            "Cc3cccc(O)c3",
            "c1ncn(C)c1",
            "CC(O)(Cc1cc2ccncc2[nH]1)C(F)(F)F",
            "Cc3ccc(F)cc3C",
            "Cc1c[nH]c2c(cccc12)",
            "CC1CC(C)(C)Nc2c(C)cc(cc12)",
            "C(CN)C(F)(F)",
            "CCOc1ccc(cc1)",
            "CC(O)(Cc1cc2ncncc2[nH]1)C(F)(F)F"
    };

    private final static String[] SmartsTertiary_2_1_1 = {
            "O=C1NC=C(Cc2ccccc2)C(=O)N1",
//            "n1c(=O)[nH]c(c(Cc2ccccc2)c1=O)", this has been changed with the prev one, structure can not be aromatic
            "N1CCC2=CC(=O)CCC2(Cc2ccccc2)C1"
    };

    private final static String[] SmartsTertiary_3_1_1 = {
            "C(C)C(F)",
            "C(O)C1CCCC2=CC(=O)CCC12C",
            "SCCO",
    };
    private final static String[] SmartsTertiary_3_1_1_1 = {
            "SC(CCO)C(=O)",
            "C(C)C(F)(F)F"
    };

    private final static String[] SmartsTertiary_3_2_1 = {
            "CC12CCC(CC1CCc1cc(O)ccc21)"
    };

    private final static String[] SmartsTertiary_3_2_1_1 = {
            "OCc3ccccc3",
            "C#CC",
            "CCC12CCC(CC1CCc1cc(O)ccc21)",
    };

    private final static String[] SmartsTertiary_3_3_1 = {
            "c1ccc(OC)cc1"
    };

    public Short getMatches() {
        return Matches;
    }

    public SAGlucocorticoidReceptor() throws InitFailureException {

        Matches = NOT_PREDICTED;

        try {
            SA_SmartsPrimary_I = new SmartsPattern[SmartsPrimaryI.length];
            SA_SmartsPrimary_II = new SmartsPattern[SmartsPrimaryII.length];
            SA_SmartsPrimary_III = new SmartsPattern[SmartsPrimaryIII.length];
            SA_SmartsSecondary_1_1 = new SmartsPattern[SmartsSecondary_1_1.length];
            SA_SmartsSecondary_2_1 = new SmartsPattern[SmartsSecondary_2_1.length];
            SA_SmartsSecondary_3_1 = new SmartsPattern[SmartsSecondary_3_1.length];
            SA_SmartsTertiary_1_1_1 = new SmartsPattern[SmartsTertiary_1_1_1.length];
            SA_SmartsTertiary_1_1_2 = new SmartsPattern[SmartsTertiary_1_1_2.length];
            SA_SmartsTertiary_1_1_3 = new SmartsPattern[SmartsTertiary_1_1_3.length];
            SA_SmartsTertiary_1_1_4 = new SmartsPattern[SmartsTertiary_1_1_4.length];
            SA_SmartsTertiary_2_1_1 = new SmartsPattern[SmartsTertiary_2_1_1.length];
            SA_SmartsTertiary_3_1_1 = new SmartsPattern[SmartsTertiary_3_1_1.length];
            SA_SmartsTertiary_3_1_1_1 = new SmartsPattern[SmartsTertiary_3_1_1_1.length];
            SA_SmartsTertiary_3_2_1 = new SmartsPattern[SmartsTertiary_3_2_1.length];
            SA_SmartsTertiary_3_2_1_1 = new SmartsPattern[SmartsTertiary_3_2_1_1.length];
            SA_SmartsTertiary_3_3_1 = new SmartsPattern[SmartsTertiary_3_3_1.length];

            // init smarts pattern

            int idx = 0;
            for (String smarts : SmartsPrimaryI) {
                SA_SmartsPrimary_I[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsPrimary_I[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsPrimaryII) {
                SA_SmartsPrimary_II[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsPrimary_II[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsPrimaryIII) {
                SA_SmartsPrimary_III[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsPrimary_III[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsSecondary_1_1) {
                SA_SmartsSecondary_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsSecondary_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsSecondary_2_1) {
                SA_SmartsSecondary_2_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsSecondary_2_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsSecondary_3_1) {
                SA_SmartsSecondary_3_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsSecondary_3_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_1_1_1) {
                SA_SmartsTertiary_1_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_1_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_1_1_2) {
                SA_SmartsTertiary_1_1_2[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_1_1_2[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_1_1_3) {
                SA_SmartsTertiary_1_1_3[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_1_1_3[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_1_1_4) {
                SA_SmartsTertiary_1_1_4[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_1_1_4[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_2_1_1) {
                SA_SmartsTertiary_2_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_2_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_3_1_1) {
                SA_SmartsTertiary_3_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_3_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_3_1_1_1) {
                SA_SmartsTertiary_3_1_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_3_1_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_3_2_1) {
                SA_SmartsTertiary_3_2_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_3_2_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_3_2_1_1) {
                SA_SmartsTertiary_3_2_1_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_3_2_1_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

            idx = 0;
            for (String smarts : SmartsTertiary_3_3_1) {
                SA_SmartsTertiary_3_3_1[idx] = SmartsPattern.create(smarts, DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_SmartsTertiary_3_3_1[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + smarts);
                idx++;
            }

        } catch (Exception ex){
            throw new InitFailureException("Unable to initialize SMARTS - " + ex.getMessage());
        }
    }

    public void Match(InsilicoMolecule mol) {

        try {

            boolean found_primary = false;

            // check primary fragments
            for(int i = 0; i < SA_SmartsPrimary_I.length; i++){
                if(SA_SmartsPrimary_I[i].matches(mol.GetStructure())) {
                    found_primary = true;

                    // check 1-1
                    boolean found_secondary = false;
                    for(int j = 0; j < SA_SmartsSecondary_1_1.length; j++){

                        if(SA_SmartsSecondary_1_1[j].matches(mol.GetStructure())){
                            found_secondary = true;

                            // check for workflow I
                            boolean found_tertiary = false;
                            for(int k = 0; k < SA_SmartsTertiary_1_1_1.length; k++){
                                if(SA_SmartsTertiary_1_1_1[k].matches(mol.GetStructure())){
                                    found_tertiary = true;
                                    Matches = ACTIVE_AGONIST;
                                }
                                if(found_tertiary)
                                    break;
                            }

                            if(!found_tertiary){
                                for(int k = 0; k < SA_SmartsTertiary_1_1_2.length; k++){
                                    if(SA_SmartsTertiary_1_1_2[k].matches(mol.GetStructure())){
                                        found_tertiary = true;
                                        Matches = ACTIVE_ANTAGONIST;
                                    }
                                    if(found_tertiary)
                                        break;
                                }
                            }

                            if(!found_tertiary){
                                for(int k = 0; k < SA_SmartsTertiary_1_1_3.length; k++){
                                    if(SA_SmartsTertiary_1_1_3[k].matches(mol.GetStructure())){
                                        found_tertiary = true;
                                        Matches = ACTIVE_AGONIST;
                                    }
                                    if(found_tertiary)
                                        break;
                                }
                            }

                            if(!found_tertiary){
                                for(int k = 0; k < SA_SmartsTertiary_1_1_4.length; k++){
                                    if(SA_SmartsTertiary_1_1_4[k].matches(mol.GetStructure())){
                                        found_tertiary = true;
                                        Matches = ACTIVE_A_ANTA;
                                    }
                                    if(found_tertiary)
                                        break;
                                }
                            }

                            if(!found_tertiary){
                                Matches = ACTIVE_ANTAGONIST;
                                break;
                            }
                        }
                    }

                    // no matches found in secondary 1-1 -> molecule is INACTIVE
                    if(!found_secondary){
                        Matches = INACTIVE;
                        break;
                    }
                }
                if(found_primary)
                    break;
            }

            if(!found_primary){
                for(int i = 0; i < SA_SmartsPrimary_II.length; i++){

                    if(SA_SmartsPrimary_II[i].matches(mol.GetStructure())){
                        found_primary = true;

                        // check 2-1
                        boolean found_secondary = false;
                        for(int j = 0; j < SA_SmartsSecondary_2_1.length; j++){

                            if(SA_SmartsSecondary_2_1[j].matches(mol.GetStructure())){
                                found_secondary = true;

                                // check for workflow II
                                boolean found_tertiary = false;

                                // search for tertiary 2-1-1 fragments
                                for(int k = 0; k < SA_SmartsTertiary_2_1_1.length; k++){

                                    if(SA_SmartsTertiary_2_1_1[k].matches(mol.GetStructure())){
                                        found_tertiary = true;
                                        Matches = ACTIVE_ANTAGONIST;
                                        break;
                                    }

                                }
                                if(!found_tertiary){
                                    Matches = ACTIVE_A_ANTA;
                                }
                            }

                        }

                        // no matches found in secondary 2-1 -> molecule is INACTIVE
                        if(!found_secondary) {
                            Matches = INACTIVE;
                            break;
                        }

                    }
                }
            }

            if(!found_primary){
                for(int i = 0; i < SA_SmartsPrimary_III.length; i++){
                    if(SA_SmartsPrimary_III[i].matches(mol.GetStructure())){
                        found_primary = true;
                        // check 3-1

                        boolean found_secondary = false;
                        for(int j = 0; j < SA_SmartsSecondary_3_1.length; j++){

                            if(SA_SmartsSecondary_3_1[j].matches(mol.GetStructure())){
                                found_secondary = true;

                                boolean found_tertiary_3_1_1 = false;

                                // check for workflow III - Search for tertiary 3-1-1
                                for(int k = 0; k < SA_SmartsTertiary_3_1_1.length; k++) {

                                    if(SA_SmartsTertiary_3_1_1[k].matches(mol.GetStructure())) {

                                        found_tertiary_3_1_1 = true;


                                        // search for tertiary 3_1_1_1
                                        boolean found_tertiary_3_1_1_1 = false;
                                        for(int y = 0; y < SA_SmartsTertiary_3_1_1_1.length; y++) {

                                            if(SA_SmartsTertiary_3_1_1_1[y].matches(mol.GetStructure())){
                                                found_tertiary_3_1_1_1 = true;
                                                Matches = ACTIVE_AGONIST;
                                            }
                                        }

                                        if(found_tertiary_3_1_1_1)
                                            break;

                                        if(!found_tertiary_3_1_1_1)
                                            Matches = ACTIVE_A_ANTA;
                                    }

                                    if(found_tertiary_3_1_1)
                                        break;

                                }


                                // search for tertiary 3_2_1
                                boolean found_tertiary_3_2_1 = false;
                                if(!found_tertiary_3_1_1){
                                    for(int k = 0; k < SA_SmartsTertiary_3_2_1.length; k++){
                                        if(SA_SmartsTertiary_3_2_1[k].matches(mol.GetStructure())){

                                            found_tertiary_3_2_1 = true;

                                            boolean found_tertiary_3_2_1_1 = false;

                                            // check for tertiary 3-2-1-1
                                            for(int y = 0; y < SA_SmartsTertiary_3_2_1_1.length; y++){
                                                if(SA_SmartsTertiary_3_2_1_1[y].matches(mol.GetStructure())){
                                                    found_tertiary_3_2_1_1 = true;
                                                    Matches = ACTIVE_A_ANTA;
                                                }
                                            }

                                            if(found_tertiary_3_2_1_1)
                                                break;

                                            if(!found_tertiary_3_2_1_1)
                                                Matches = ACTIVE_ANTAGONIST;

                                        }
                                    }
                                }


                                boolean found_tertiary_3_3_1 = false;
                                // search for tertiary 3_3_1
                                if ( (!found_tertiary_3_1_1) && (!found_tertiary_3_2_1) ) {
                                    for (int k = 0; k < SA_SmartsTertiary_3_3_1.length; k++){
                                        if(SA_SmartsTertiary_3_1_1[k].matches(mol.GetStructure())){
                                            found_tertiary_3_3_1 = true;
                                            Matches = ACTIVE_ANTAGONIST;
                                        }
                                        if(found_tertiary_3_3_1)
                                            break;
                                    }
                                }

                                if ( (!found_tertiary_3_1_1) && (!found_tertiary_3_2_1) && (!found_tertiary_3_3_1) )
                                    Matches = ACTIVE_AGONIST;
                            }
                        }

                        // no matches found in secondary 3-1 -> molecule is INACTIVE
                        if(!found_secondary){
                            Matches = INACTIVE;
                            break;
                        }
                    }
                }
            }

            if(!found_primary)
                Matches = INACTIVE;

        } catch (Exception ex){
            log.warn(ex.getMessage());
        }
    }


}
