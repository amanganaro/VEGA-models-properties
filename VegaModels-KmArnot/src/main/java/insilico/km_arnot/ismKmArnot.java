package insilico.km_arnot;

import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;

import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;

import insilico.core.tools.utils.ModelUtilities;
import insilico.km_arnot.descriptors.KmFactor;
import insilico.km_arnot.descriptors.MeylanLogP;
import insilico.km_arnot.descriptors.ad.ADIndexADILogP;
import insilico.km_arnot.descriptors.weights.Mass;
import insilico.meylanlogp.ismLogPMeylan;
import org.openscience.cdk.interfaces.IAtom;

import java.util.ArrayList;

/**
 *
 * @author User
 */
public class ismKmArnot extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_km_arnot.xml";
    
    private final ismLogPMeylan LogPModel;
    
    private InsilicoModelOutput LogPResult;
    private double LogPValueForModel;

    
    public ismKmArnot() throws InitFailureException {
        super(ModelData);
        
        // Builds model objects
        LogPModel = new ismLogPMeylan();
        
        // Define no. of descriptors
        this.DescriptorsSize = 3;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "MW";
        this.DescriptorsNames[1] = "LogP";
        this.DescriptorsNames[2] = "Factors";
        
        // Defines results
        this.ResultsSize = 5;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted log kM/Half-Life (days)";
        this.ResultsName[1] = "Predicted kM/Half-Life (days)";
        this.ResultsName[2] = "Predicted LogP (Meylan/Kowwin)";
        this.ResultsName[3] = "Predicted LogP reliability";
        this.ResultsName[4] = "MW";
        
        // Define AD items
        this.ADItemsName = new String[4];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexMaxError().GetIndexName();
        
    }
    
    
    @Override
    public ArrayList<DescriptorBlock> GetRequiredDescriptorBlocks() {
        
        ArrayList<DescriptorBlock> blocks = new ArrayList<>();
        DescriptorBlock desc;

        // LogP block required, it is used in the Meylan LogP Model
        desc = new MeylanLogP();
        blocks.add(desc);

        // kM factor
        desc = new KmFactor();
        blocks.add(desc);

        return blocks;
    }    
    
    
    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {
            Descriptors = new double[DescriptorsSize];
            
            // Calculate the whole Meylan LogP model
            LogPResult = LogPModel.Execute(CurMolecule, null, false);
            if (LogPResult.HasExperimental())
                LogPValueForModel = LogPResult.getExperimental();
            else
                LogPValueForModel = LogPResult.getMainResultValue();
            
            // Calculates correct MW
            double[] wMass = Mass.getWeights(CurMolecule.GetStructure());
            int H = 0;
            for (int i=0; i<CurMolecule.GetStructure().getAtomCount(); i++) { 
                IAtom CurAt = CurMolecule.GetStructure().getAtom(i);
                try {
                    H += CurAt.getImplicitHydrogenCount();
                } catch (Exception e) { }
            }
            
            double MW = H * 1.010;
            for (double mass : wMass) {
                if (mass == -999) {
                    MW = -999;
                    break;
                }
                MW += mass * 12.0107;
            }
            
            Descriptors = new double[DescriptorsSize];
            Descriptors[0] = MW;
            Descriptors[1] = LogPValueForModel;
            Descriptors[2] = DescEngine.GetDescriptorBlock(KmFactor.class).GetByName("kM").getValue();
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {

        if (LogPResult.getStatus() != InsilicoModelOutput.OUTPUT_OK)
            return MODEL_ERROR;
        
        // Arnot formula:
        // HL = 0.30734215 * LogP - 0.0025643319 * MW - 1.53706847 + FragCoeff;        
        double HL = 0.30734215 * Descriptors[1] - 0.0025643319 * Descriptors[0] - 1.53706847 + Descriptors[2];        
        
        CurOutput.setMainResultValue(HL);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(HL)); // log km/HL
        Res[1] = String.valueOf(Format_2D.format(Math.pow(10, HL) )); // km/HL
        Res[2] = String.valueOf(Format_2D.format(LogPValueForModel)); // LogP
        String Rel = "n.a.";
        if (LogPResult.HasExperimental())
            Rel = "Experimental";
        else
            switch (LogPResult.getADI().GetAssessmentClass()) {
                case ADIndex.INDEX_LOW:
                    Rel = "Low";
                    break;
                case ADIndex.INDEX_MEDIUM:
                    Rel = "Moderate";
                    break;
                case ADIndex.INDEX_HIGH:
                    Rel = "Good";
                    break;
            }
        Res[3] = Rel; // Reliability of logp prediction
        Res[4] = String.valueOf(Format_2D.format(Descriptors[0])); // MW
        
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
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.0, 0.5);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.0, 0.5);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.0, 0.5);
        } catch (Throwable e) {
            return InsilicoModel.AD_ERROR;
        }
        
        // Sets final AD index
        ADIndexADI ADI = new ADIndexADILogP(adq.getIndexADI(), CurOutput.getADIndex(ADIndexAccuracy.class),
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
