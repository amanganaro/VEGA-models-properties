package insilico.verhaar_toxtree;

import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import lombok.extern.log4j.Log4j;
import shadedTTv31.org.openscience.cdk.DefaultChemObjectBuilder;
import shadedTTv31.org.openscience.cdk.exception.InvalidSmilesException;
import shadedTTv31.org.openscience.cdk.interfaces.IAtomContainer;
import shadedTTv31.org.openscience.cdk.smiles.SmilesParser;
import toxTree.core.IDecisionMethod;
import toxTree.core.IDecisionResult;
import toxTree.exceptions.DecisionMethodException;
import toxTree.exceptions.DecisionResultException;
import toxTree.exceptions.MolAnalyseException;
import toxTree.query.MolAnalyser;
import verhaar.VerhaarScheme;


/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class ismVerhaarToxtree extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_verhaar_toxtree.xml";
        
    private final IDecisionMethod ttVerhaar;
    private final SmilesParser SP;
    
    
    public ismVerhaarToxtree() 
            throws InitFailureException {
        super(ModelData);
        
        // Init TT module and its smiles parser
        try {
            ttVerhaar = new VerhaarScheme();
            SP = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        } catch (DecisionMethodException e) {
            throw new InitFailureException("Unable to init toxtree module - " + e.getMessage());
        }

        // Define no. of descriptors
        this.DescriptorsSize = 0;
        this.DescriptorsNames = new String[DescriptorsSize];
        
        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Verhaar class";
        
        // Define AD items
        this.ADItemsName = new String[0];
        
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
        
        int MainResult = -1;
        String MainResultText = "N.A.";
        
        try {
            
            // Create result object
            IDecisionResult ttResult = ttVerhaar.createDecisionResult();

            // Parse the input molecule
            IAtomContainer ac = SP.parseSmiles(this.CurMolecule.GetSMILES());

            // Run the tree
            MolAnalyser.analyse(ac);
            ttResult.classify(ac);
            ttResult.assignResult(ac);

            // Assumes there only one assigned final category available
            MainResult = ttResult.getAssignedCategories().get(0).getID();
            MainResultText = ttResult.getAssignedCategories().get(0).getName();
            
        } catch (InvalidSmilesException | DecisionResultException | MolAnalyseException e) {
            //
        }

        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        Res[0] = MainResultText;
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {
        
        return InsilicoModel.AD_ERROR;
        
    }
        

    @Override
    protected void CalculateAssessment() {
        
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);
        
    }
    
}
