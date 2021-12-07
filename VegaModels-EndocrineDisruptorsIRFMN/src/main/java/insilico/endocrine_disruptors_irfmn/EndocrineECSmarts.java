package insilico.endocrine_disruptors_irfmn;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;

public class EndocrineECSmarts {

    private SmartsPattern[] SA;
    private SmartsPattern SA_add_3;

    final static String[][] SMARTS = {
            {"[$(C[Cl]);!$(C(Cl)(Cl)Cl)]", "Polychlorinated (aliphatic)"}, // 0
            {"[$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)]-[$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)]", "PCBs"},
            {"[$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)][O,o][$(c1cccc([Cl,Br])c1[Cl,Br]),$(c1ccc([Cl,Br])cc1[Cl,Br]),$(c1cc([Cl,Br])ccc1[Cl,Br]),$(c1c([Cl,Br])cccc1[Cl,Br]),$(c1ccc([Cl,Br])c([Cl,Br])c1),$(c1cc([Cl,Br])cc([Cl,Br])c1),$(c1ccc([Cl,Br])c([Cl,Br])c1)]", "PCDDs"},
            {"c1ccccc1Cc2ccccc2", "Polyhalogenated aromatic compounds"}, // 3
            {"S=C(N)(S)", "S=C(N)(S) group"},
            {"O(C)C(=O)c1ccccc1C(OC)=O", "Phthalates"},
            {"[Sn]([C,c])([C,c])[C,c]", "Tin compounds"},
            {"Oc1ccccc1CCCCCC", "Phenol with aliphatic chain"},
            {"Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C", "Bisphenol A like compounds"}
    };

    public EndocrineECSmarts() throws InitFailureException {
        try {
            SA = new SmartsPattern[SMARTS.length];
            int idx = 0;
            for (String[] s : SMARTS) {
                SA[idx] = SmartsPattern.create(s[0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                if (SA[idx] == null)
                    throw new Exception("unable to parse SMARTS: " + s[0]);
                idx++;
            }
            SA_add_3 = SmartsPattern.create("[Cl,Br]", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS");
        }

    }

    public String Match(InsilicoMolecule mol) throws GenericFailureException {

        try {

            for (int i=0; i<SA.length; i++)
                if (SA[i].matches(mol.GetStructure())) {
                    if (i==0) {
                        int count = SA[i].matchAll(mol.GetStructure()).countUnique();
                        if (count > 5)
                            return(SMARTS[i][1]);
                    } else if (i==3) {
                        if (SA_add_3.matches(mol.GetStructure())) {
                            if (SA_add_3.matchAll(mol.GetStructure()).countUnique() >= 4)
                                return(SMARTS[i][1]);
                        }
                    } else
                        return(SMARTS[i][1]);
                }

        } catch (Exception ex) {
            throw new GenericFailureException(ex.getMessage());
        }

        return null;
    }


}
