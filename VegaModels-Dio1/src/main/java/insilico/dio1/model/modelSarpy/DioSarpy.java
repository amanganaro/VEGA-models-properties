package insilico.dio1.model.modelSarpy;

import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;

public class DioSarpy {

    private final static String[][] RULES = {
            {"CCCCCCCC[N+](C)(C)C","active","inf"},
            {"CCOCCOC","inactive","inf"},
            {"C(N)C(=O)N","inactive","inf"},
            {"O=C(O)c1ccccc1C","inactive","inf"},
            {"CC1(CCCC1)C(C)","inactive","inf"},
            {"Oc1ccc(cc1)CCCC","active","8.15"},
            {"c4ccc(cc4)c4ccccc4","inactive","inf"},
            {"P(OCC)OC","inactive","inf"},
            {"P(=S)(OC)OC","active","2.83"},
            {"Oc1ccccc1OC","inactive","inf"},
            {"N(CNC)C","inactive","inf"},
            {"C(OCCc1ccccc1)","inactive","inf"},
            {"CCCCCCCC=CC=CC","active","8.49"},
            {"C(C)CC(C)(C)C","inactive","1.91"},
            {"n1c(nc(nc1))","inactive","inf"},
            {"COc1ccc(cc1)Cl","inactive","inf"},
            {"CCNCCN","inactive","inf"},
            {"C(OCCCCCC)CCCCC","inactive","inf"},
            {"ON=C","inactive","inf"},
            {"NN(C)","inactive","inf"},
            {"Cc1cccc(Oc2ccccc2)c1","inactive","inf"},
            {"Cc1ncccc1","inactive","inf"},
            {"C1CC2C=CC1C2","inactive","inf"},
            {"C(N(C)C)S","inactive","inf"},
            {"N(O)C","inactive","inf"},
            {"Nc1cccc2ccccc12","inactive","inf"},
            {"Cc1ccc(cc1)C(C)(C)C","inactive","inf"},
            {"O=C(c1cccc(c1))N(CC)","inactive","inf"},
            {"Cn1ncnc1","inactive","inf"},
            {"O=C(Nc1cccc(c1))N(C)","inactive","inf"},
            {"O=c1ncccn1","inactive","inf"},
            {"O=c1c(O)c(oc2cc(O)cc(O)c12)c1ccc(O)cc1","active","inf"},
            {"O=C(N(c1c(cccc1CC)C)C)CCl","active","inf"},
            {"C(c1ccc(cc1))C1CCCCC1","active","inf"},
            {"c1cc(c(c(c1NC(CC))))C","inactive","inf"},
            {"O=C(OCc1ccccc1)","inactive","inf"},
            {"O=C(O)CCCCC(=O)","inactive","inf"},
            {"C(N)c1cnccc1","inactive","inf"},
            {"OC(OCCC)","inactive","inf"},
            {"C1COCC1","inactive","inf"},
            {"CC(C)CCCC(O)(C)C","inactive","inf"},
            {"O=C(OCCCC)C(O)","inactive","inf"},
            {"c1cc(nnc1)","inactive","inf"},
            {"[Si](O)(O)","inactive","inf"},
            {"Nc1ccc(cc1)Cc1ccc(N)c(c1)","inactive","inf"},
            {"Oc1ccc(Oc2ccccc2)cc1","inactive","inf"},
            {"O=S(=O)(O)c1ccc(cc1)C","inactive","inf"},
            {"c1cccc2cccnc12","inactive","inf"},
            {"CC=C(C)CCC=C(C)C","inactive","inf"},
            {"O=C(c1cccc(c1)C)","inactive","inf"},
            {"O=C(Oc1ccccc1)","inactive","inf"},
            {"Nc1ccc(c(c1)Cl)Cl","inactive","inf"},
            {"O(c1cc(ccc1N))C","inactive","inf"},
            {"O=CCC(O)CO","inactive","inf"},
            {"O=c1oc(ccc1)","inactive","inf"},
            {"O(CCN(C)C)C","inactive","inf"},
            {"C(OC)NCCC","inactive","inf"},
            {"CCBr","inactive","inf"},
            {"O=Cc1ccccc1Cc1ccccc1","inactive","inf"},
            {"O=C(OC)c1ccc(cc1)C(=O)OC","inactive","inf"},
            {"C(C)C(c1ccccc1)CC(N(C))","inactive","inf"},
            {"O=C(OC(C)(C)CCC=C(C)C)C","inactive","inf"},
            {"c1cc(cc2ccc(cc12)C)C","inactive","inf"},
            {"O=P(Oc1ccccc1)O","inactive","inf"},
            {"C(N)Cc1ccc(O)cc1","inactive","inf"},
            {"c1ccc(O)c(c1)[N+](=O)","inactive","inf"},
            {"Nc1ccc(c(N)c1)C","inactive","inf"},
            {"Nc1ccc(cc1)CC","inactive","inf"},
            {"n1nc2ccccc2n1","inactive","inf"},
            {"Cn3cncc3","inactive","inf"},
            {"C(O)C(Cl)(Cl)","inactive","inf"},
    };

    private SmartsPattern[] SAS;

    public DioSarpy() {
        SAS = new SmartsPattern[RULES.length];
        for (int i=0; i<RULES.length; i++) {
            SmartsPattern sp = SmartsPattern.create(RULES[i][0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            SAS[i] = sp;
        }
    }

    public int match(InsilicoMolecule mol) throws Exception {

        for (int i=0; i<RULES.length; i++) {
            if (SAS[i].matches(mol.GetStructure())) {
                if (RULES[i][1].equalsIgnoreCase("active"))
                    return 1;
                else if (RULES[i][1].equalsIgnoreCase("inactive"))
                    return 0;
                else
                    throw new Exception ("wrong field active/inactive in SMARTS list");
            }
        }
        return -1;
    }

}
