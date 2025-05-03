package insilico.ld50quail.alerts;

import insilico.core.alerts.*;
import insilico.core.constant.InsilicoConstants;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.tools.Depiction;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.ArrayList;

/**
 *
 * @author User
 */
@Slf4j
public class SALD50Quail extends AlertBlockFromSMARTS implements iAlertBlock {

    private Pattern[] SA;

    private final static String[][]SMARTS_L = {
        // all inf
        { "P(OCC)","ACTIVE" },
        { "O=CC(=COC)c1ccccc1","INACTIVE" },
        { "c2ccc(cc2)Br","ACTIVE" },
        { "P(=O)","ACTIVE" },
        { "C(c1ccc(cc1)Cl)C(C)C","ACTIVE" },
        { "C(CNC(c1ccc(cc1)))C","INACTIVE" },
        { "c1nc(OC)cc(n1)OC","INACTIVE" },
        { "Fc1ccc(cc1)","INACTIVE" },
        { "NN(C)C(C)","INACTIVE" },
        { "c2ccc(cc2)c2ccccc2","ACTIVE" },
        { "O=C(O)COc1ccc(cc1)Cl","ACTIVE" },
        { "Oc1cccc(c1)Cl","ACTIVE" },
        { "C(N(CCOC))","ACTIVE" },
        { "C(C(Cl))Cl","ACTIVE" },
        { "O=C(Oc1cc(c(c(c1))))NC","ACTIVE" },
        { "O=[N+]([O-])c1ccc(O)cc1","ACTIVE" },
        { "Nc1ccc(c(c1)Cl)","ACTIVE" },
        { "CSP(=S)(OC)OC","ACTIVE" },
        { "SCCC","ACTIVE" },
        { "CC#C","ACTIVE"}
    };

    private final static String[][] SMARTS_H = {
        // all inf
        {"O=CCC", "INACTIVE"},
        {"P(=O)(OC)", "ACTIVE"},
        {"C(C)COC", "INACTIVE"},
        {"NCN", "INACTIVE"},
        {"OP(=S)(OCC)OCC", "ACTIVE"},
        {"Oc2ccc(cc2Cl)", "INACTIVE"},
        {"OP(OC)OC", "ACTIVE"},
        {"c1cncnc1", "INACTIVE"},
        {"C(c1ccccc1)CC(O)", "INACTIVE"},
        {"N(CCC)C", "INACTIVE"},
        {"O=C(Oc1cc(c(c(c1))))NC", "ACTIVE"},
        {"O=CCS", "ACTIVE"},
        {"Cc1ccccc1Cl", "INACTIVE"},
        {"c1cc(ccc1C)C", "INACTIVE"},
        {"n1c(nc(nc1))NC", "INACTIVE"},
    };


    public SALD50Quail() throws InitFailureException {
        super(1, "Rules for LD50 in quail (IRFMN)");
    }
    
    
    @Override
    protected void BuildSAList() throws InitFailureException {

        int idx = 0;

        for (int i=0; i<SMARTS_L.length; i++) {
            Alert curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (idx+1)));
            curSA.setName("LD50 (quail) Low alert no. " + (i+1));
            curSA.setDescription("Structural alert for LD50 (quail) Low toxicity - " + SMARTS_L[i][1] + " defined by the SMARTS: " + SMARTS_L[i][0]);
            curSA.setBoolProperty("LD50QUAIL_ACTIVE", SMARTS_L[i][1].equalsIgnoreCase("ACTIVE") );
            curSA.setBoolProperty("LD50QUAIL_L", true);
            Alerts.add(curSA);
            idx++;
        }

        for (int i=0; i<SMARTS_H.length; i++) {
            Alert curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (idx+1)));
            curSA.setName("LD50 (quail) High alert no. " + (i+1));
            curSA.setDescription("Structural alert for LD50 (quail) High toxicity - " + SMARTS_H[i][1] + " defined by the SMARTS: " + SMARTS_H[i][0]);
            curSA.setBoolProperty("LD50QUAIL_ACTIVE", SMARTS_H[i][1].equalsIgnoreCase("ACTIVE") );
            curSA.setBoolProperty("LD50QUAIL_H", true);
            Alerts.add(curSA);
            idx++;
        }

    }
    
    
    @Override
    protected void InitSMARTS() throws InitFailureException {

        try {
            int nFragments = SMARTS_L.length + SMARTS_H.length;
            SA = new Pattern[nFragments];
            
            int idx = 0;
            for (String s[] : SMARTS_L) {
                SA[idx] = SmartsPattern.create(s[0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                idx++;
            }
            for (String s[] : SMARTS_H) {
                SA[idx] = SmartsPattern.create(s[0], DefaultChemObjectBuilder.getInstance()).setPrepare(false);
                idx++;
            }
        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS");
        }    
    }

    
    @Override
    protected AlertList CalculateSAMatches() throws GenericFailureException {
        AlertList Res = new AlertList();
        
        try {

            int nFragments = SMARTS_L.length + SMARTS_H.length;
            
            for (int i=0; i<nFragments; i++) 
                if ((SA[i].matches(CurMol.GetStructure())))
                    Res.add((Alert)Alerts.get(i).clone());
            
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            return null;
        }
        
        return Res; 
    }

}