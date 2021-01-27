package insilico.skin_permeation_potts;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.tools.utils.ModelUtilities;
import insilico.logp_alogp.ismLogPALogP;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import static insilico.skin_permeation_potts.weights.Mass.GetMass;


/**
 *
 * @author User
 */
public class ismSkinPermeationPotts extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_skin_perm_potts.xml";

    private double ALogP;
    private double MW;
    private final InsilicoModel mAlogP;
    
    
    public ismSkinPermeationPotts() 
            throws InitFailureException {
        super(ModelData);
        
        // Init logP model
        mAlogP = new ismLogPALogP();
        
        // Define no. of descriptors
        this.DescriptorsSize = 2;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "MW";
        this.DescriptorsNames[1] = "ALogP";
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted skin permeation logKp [log(cm/h)]";
        this.ResultsName[1] = "Molecular Weight";
        this.ResultsName[2] = "AlogP";
        
        // Define AD items
        this.ADItemsName = new String[6];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        this.ADItemsName[4] = new ADIndexRange().GetIndexName();
        this.ADItemsName[5] = new ADIndexACF().GetIndexName();

    }
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            
            Descriptors = new double[DescriptorsSize];
            
            double MASS_C = 12.0107;
            double MASS_H = 1.00784;
            
            // MW is calculated here to be more precise (bug to be fixed in the core
            // libraries, because we have the weights for relative atomic mass)
            double MW = 0;
            IAtomContainer curMol = CurMolecule.GetStructure();
            int nSK = curMol.getAtomCount();
            for (int i=0; i<nSK; i++) {
                IAtom CurAt = curMol.getAtom(i);
                double m = GetMass(CurAt.getSymbol());
                if (m == Descriptor.MISSING_VALUE)
                    throw new InvalidMoleculeException();
                MW += (m * MASS_C);

                // Hydrogens
                int H = 0;
                try {
                    H = CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
                MW += (H * MASS_H);
            }            
            

            // logP
            InsilicoModelOutput o = mAlogP.Execute(CurMolecule);
            double logp;
            if (o.HasExperimental())
                logp = o.getExperimental();
            else 
                if (o.getStatus() >= InsilicoModelOutput.OUTPUT_OK)
                    logp = o.getMainResultValue();
                else
                    throw new Exception("Error in AlogP model");
            
            
            Descriptors[0] = MW;
            Descriptors[1] = logp;

            this.MW = MW;
            this.ALogP = logp;
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {
        
        // Potts and Guy formula
        double formula = 0.71 * Descriptors[1] - 0.0061 * Descriptors[0] - 2.7;
        
        CurOutput.setMainResultValue(formula);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_3D.format(formula)); // log(cm/h)
        Res[1] = Format_2D.format(MW);
        Res[2] = Format_3D.format(ALogP);
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQuantitative adq = new ADCheckIndicesQuantitative(TS);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.85, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.2, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.2, 0.6);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.2, 0.6);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
        // Sets Range check
        ADCheckDescriptorRange adrc = new ADCheckDescriptorRange();
        if (!adrc.Calculate(TS, Descriptors, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets ACF check
        ADCheckACF adacf = new ADCheckACF(TS);
        if (!adacf.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;
        
        // Sets final AD index
        double acfContribution = CurOutput.getADIndex(ADIndexACF.class).GetIndexValue();
        double rcContribution = CurOutput.getADIndex(ADIndexRange.class).GetIndexValue();
        double ADIValue = adq.getIndexADI() * acfContribution * rcContribution;

        ADIndexADIAggregate ADI = new ADIndexADIAggregate(0.85, 0.7, 1, 0.85, 0.7);
        ADI.SetValue(ADIValue, CurOutput.getADIndex(ADIndexAccuracy.class),
                CurOutput.getADIndex(ADIndexConcordance.class),
                CurOutput.getADIndex(ADIndexMaxError.class));
        CurOutput.setADI(ADI);
        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
       
        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Always gray light - no threshold for this model
        CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GRAY);

    }
}
