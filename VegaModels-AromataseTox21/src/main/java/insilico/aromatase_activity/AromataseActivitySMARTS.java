package insilico.aromatase_activity;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.ArrayList;

public class AromataseActivitySMARTS {

    private final static String ACTIVE = "Active";
    private final static String INACTIVE = "Inactive";
    private final static String ACTIVE_AGONIST = "Active agonist";
    private final static String ACTIVE_ANTAGONIST = "Active antagonist";

    private final static String[][] SMARTS_ACTIVITY = {
            {"n1c2ccc(OC)cc2sc1N",ACTIVE,"SA1"},
            {"c1ccc(cc1)c1nc(N)sc1",ACTIVE,"SA2"},
            {"n1c(N)sc2cccc(c12)",ACTIVE,"SA3"},
            {"CCCCN",INACTIVE,"SA4"},
            {"O=S(=O)(N)c1ccc(N)cc1",INACTIVE,"SA5"},
            {"n1ccsc1C",INACTIVE,"SA6"},
            {"O=C(O)CC",INACTIVE,"SA7"},
            {"N=C(N)N",INACTIVE,"SA8"},
            {"n1oc(cc1)",INACTIVE,"SA9"},
            {"n1csc(c1C)",INACTIVE,"SA10"},
            {"O=c1oc2cc(ccc2[nH]1)",INACTIVE,"SA11"},
            {"c1c[n+](cn1CCCCCCCC)C",ACTIVE,"SA12"},
            {"c1c(cccc1Cl)Cl",ACTIVE,"SA13"},
            {"c1ccnn1",INACTIVE,"SA14"},
            {"c1cnc(n1C)C",INACTIVE,"SA15"},
            {"O=c1c2c(ncn2C)n(c(=O)n1)C",INACTIVE,"SA16"},
            {"c1nc[nH]c1C",INACTIVE,"SA17"},
            {"O=C(O)C",INACTIVE,"SA18"},
            {"S(=O)c1ccccc1",INACTIVE,"SA19"},
            {"O=c1ccnc([nH]1)",INACTIVE,"SA20"},
            {"O=C(NC)C",INACTIVE,"SA21"},
            {"O=C(OC)c1cncn1",INACTIVE,"SA22"},
            {"N(C)C",INACTIVE,"SA23"},
            {"OC(C)Cn1cnc2c(ncnc12)N",INACTIVE,"SA24"},
            {"C(OCC)C(C)",ACTIVE,"SA25"},
            {"N#C",ACTIVE,"SA26"},
            {"n1nncc1",INACTIVE,"SA27"},
            {"N(CC)C",INACTIVE,"SA28"},
            {"O=S(=O)(N)",INACTIVE,"SA29"}
    };

    private final static String[][] SMARTS_ACTIVITY_TYPE = {
            {"c2cscn2",ACTIVE_AGONIST,"SA1"},
            {"Clc1ccc(CC)cc1",ACTIVE_ANTAGONIST,"SA2"},
            {"n1cncn1",ACTIVE_ANTAGONIST,"SA3"},
            {"C(=O)O",ACTIVE_AGONIST,"SA4"},
            {"CCCN(CC)",ACTIVE_ANTAGONIST,"SA5"},
            {"c1cn(cn1)C(c1ccccc1)",ACTIVE_ANTAGONIST,"SA6"},
            {"c1c(cccc1)n1ccnc1",ACTIVE_ANTAGONIST,"SA7"},
            {"C(=O)N",ACTIVE_AGONIST,"SA8"},
            {"Clc1cccc(c1)",ACTIVE_ANTAGONIST,"SA9"},
            {"c1ncnc2c1ncn2",ACTIVE_AGONIST,"SA10"}
    };

    private SmartsPattern[] SA_ACTIVITY;
    private SmartsPattern[] SA_ACTIVITY_TYPE;

    public ArrayList<String> Matches_Active;
    public ArrayList<String> Matches_Inactive;
    public ArrayList<String> Matches_Active_Agonist;
    public ArrayList<String> Matches_Active_Antagonist;



    public AromataseActivitySMARTS() throws InitFailureException {

        try {
            SA_ACTIVITY = new SmartsPattern[SMARTS_ACTIVITY.length];
            SA_ACTIVITY_TYPE = new SmartsPattern[SMARTS_ACTIVITY_TYPE.length];

            int idx = 0;
            for (String[] s : SMARTS_ACTIVITY) {
                SA_ACTIVITY[idx] = SmartsPattern.create(s[0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_ACTIVITY[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + s[0]);
                idx++;
            }

            idx = 0;
            for (String[] s : SMARTS_ACTIVITY_TYPE) {
                SA_ACTIVITY_TYPE[idx] = SmartsPattern.create(s[0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA_ACTIVITY_TYPE[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + s[0]);
                idx++;
            }

        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS - " + e.getMessage());
        }

    }


    public void Match(InsilicoMolecule mol) throws GenericFailureException {

        Matches_Active = new ArrayList<>();
        Matches_Inactive = new ArrayList<>();
        Matches_Active_Agonist = new ArrayList<>();
        Matches_Active_Antagonist = new ArrayList<>();

        try {

            for (int i=0; i<SMARTS_ACTIVITY.length; i++)
                if (SA_ACTIVITY[i].matches(mol.GetStructure())) {
                    if (SMARTS_ACTIVITY[i][1].equalsIgnoreCase(ACTIVE))
                        Matches_Active.add(SMARTS_ACTIVITY[i][2]);
                    if (SMARTS_ACTIVITY[i][1].equalsIgnoreCase(INACTIVE))
                        Matches_Inactive.add(SMARTS_ACTIVITY[i][2]);
                }

            for (int i=0; i<SMARTS_ACTIVITY_TYPE.length; i++)
                if (SA_ACTIVITY_TYPE[i].matches(mol.GetStructure())) {
                    if (SMARTS_ACTIVITY_TYPE[i][1].equalsIgnoreCase(ACTIVE_AGONIST))
                        Matches_Active_Agonist.add(SMARTS_ACTIVITY_TYPE[i][2]);
                    if (SMARTS_ACTIVITY_TYPE[i][1].equalsIgnoreCase(ACTIVE_ANTAGONIST))
                        Matches_Active_Antagonist.add(SMARTS_ACTIVITY_TYPE[i][2]);
                }

        } catch (Exception ex) {
            throw new GenericFailureException(ex.getMessage());
        }

    }

    public static String FormatAlertArray(ArrayList<String> SAs) {
        String s = "";
        for (String SA : SAs) {
            if (!s.isEmpty()) s += ", ";
            s += SA;
        }
        return s;
    }

}
