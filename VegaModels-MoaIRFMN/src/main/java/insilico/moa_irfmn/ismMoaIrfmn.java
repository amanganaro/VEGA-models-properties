package insilico.moa_irfmn;

import insilico.core.ad.ADCheckSA;
import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.builders.SAMoaIRFMN;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.similarity.SimilarMolecule;
import insilico.core.tools.utils.ModelUtilities;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ismMoaIrfmn extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_moa_irfmn.xml";
     
    private final String FullAlertSet;
    private final SAMoaIRFMN SAMoA;
    
    
    public ismMoaIrfmn() 
            throws InitFailureException {
        super(ModelData);
        
        SAMoA = new SAMoaIRFMN();
        
        // Set SA list
        FullAlertSet = AlertEncoding.MergeAlertIds(SAMoA.getAlerts());
        
        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 2;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted MoA";
        this.ResultsName[1] = "Other MoA identified";
        
        // Define AD items
        this.ADItemsName = new String[0];
        
    }
    
    @Override
    public ArrayList<Integer> GetRequiredAlertBlocks() {
        ArrayList<Integer> blocks = new ArrayList<>();
//        blocks.add(InsilicoConstants.SA_BLOCK_MOA_IRFMN);
        blocks.add(51);
        return blocks;
    }    
    
    
    @Override
    public AlertList GetCalculatedAlert() throws CloneNotSupportedException {
        AlertList FoundSAs = new AlertList();
        for (Alert a : CurMolecule.GetAlerts().getSAList()) 
            if (AlertEncoding.ContainsAlert(FullAlertSet, a.getId())) {
                FoundSAs.add((Alert)a.clone());
            }
        return FoundSAs;
    }

    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        int MainResult = 0;
        String MainMoA = "No MoA";
        String AdditionalMoA = "";
        
        for (Alert a : CurMolecule.GetAlerts().getSAList()) 
            if (AlertEncoding.ContainsAlert(FullAlertSet, a.getId())) {
                if (MainResult == 0) {
                    MainResult = 1;
                    MainMoA = SAMoA.MoAList.get(a.getId());
                } else {
                    AdditionalMoA = AdditionalMoA.isEmpty() ? "" : "; ";
                    AdditionalMoA += SAMoA.MoAList.get(a.getId());
                }
            }

        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        Res[0] = MainMoA;
        Res[1] = AdditionalMoA;
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        try {
            SimilarMolecule[] sim = new SimilarMolecule[0];
            ADCheckSA ADSA = new ADCheckSA(TS);
            ADSA.Calculate(CurMolecule, GetCalculatedAlert(), CurOutput, sim);
        
        } catch (Throwable e) { }
        
        return InsilicoModel.AD_ERROR;
        
    }
        

    @Override
    protected void CalculateAssessment() {
        
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        if (CurOutput.getMainResultValue() == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        
    }
    
}
