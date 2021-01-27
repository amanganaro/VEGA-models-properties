package insilico.bcf_meylan;

import insilico.core.alerts.*;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

/**
 *
 * @author amanganaro
 */
public class BCFMeylanIonicChecker extends AlertBlockFromSMARTS implements iAlertBlock {

    private boolean IsIonic; 
    private String Message;
    private String[] SMARTSFragments;
    private Pattern[] SA;
    
    
    public BCFMeylanIonicChecker() throws Exception {
        super(-1, "Meylan BCF Ionic Fragments"); // id = -1 as it is not needed (it is not a public alert block)
    }
    
    
    @Override
    protected void BuildSAList() throws InitFailureException {

        SMARTSFragments = new String[4];
        
        Alert curSA;
        int i=0;
        
        SMARTSFragments[i] = "*C(=O)[OH]";
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFIon" + (i+1));
        curSA.setDescription("the compound contains at least one carboxylic group");
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = "*S(=O)(=O)[OH]";
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFIon" + (i+1));
        curSA.setDescription("the compound contais at least one sulfonic group");
        Alerts.add(curSA);
        i++;
            
        SMARTSFragments[i] = "*S(=O)(=O)[O-]";
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFIon" + (i+1));
        curSA.setDescription("the compound contains at least one sulfonic group in salt form");
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "[N+;!$([N+](=O)[O-])]";
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFIon" + (i+1));
        curSA.setDescription("the compound contains at least one charged nitrogen");
        Alerts.add(curSA);

    }      

    
    @Override
    protected void InitSMARTS() throws InitFailureException {

        SA = new Pattern[SMARTSFragments.length];

        try {

            for (int i=0; i<SMARTSFragments.length; i++)
                SA[i] = SmartsPattern.create(SMARTSFragments[i], DefaultChemObjectBuilder.getInstance());

        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS");
        }    
    }

    
    @Override
    protected AlertList CalculateSAMatches() throws GenericFailureException {
        
        AlertList Res = new AlertList();
        IsIonic = false;
        Message = "";
        
        try {

            for (int i=0; i<SA.length; i++) {
                
                if ((SA[i].matches(CurMol.GetStructure()))) {
                    Res.add((Alert)Alerts.get(i).clone());
                    if (!Message.isEmpty())
                        Message += "; ";
                    Message += Alerts.get(i).getDescription();
                    IsIonic = true;
                }
            }
            
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            return null;
        }
        
        return Res; 
    }    

    
    public String GetMessage() {
        return Message;
    }

    public boolean IsIonic() {
        return IsIonic;
    }
}
