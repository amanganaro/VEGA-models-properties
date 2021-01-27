package insilico.carcinogenicity_sforegression;

import insilico.carcinogenicity_sforegression.descriptors.EmbeddedDescriptors;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.tools.CustomQueryMatcher;
import insilico.core.pmml.ModelANNFromPMML;
import insilico.core.tools.utils.ModelUtilities;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author User
 */
public class ismCarcinogenicitySFORegression extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_carc_sforegression.xml";
    
    private ModelANNFromPMML ANN;
    
    
    public ismCarcinogenicitySFORegression() 
            throws InitFailureException {
        super(ModelData);
        
        // Init PMML ANN
        try {  
            URL src = getClass().getResource("/data/sfo_ann.pmml");
            ANN = new ModelANNFromPMML(src.openStream(), "Log SF");
        } catch (IOException ex) {
            throw new InitFailureException("Unable to read PMML source from .jar file");
        }
        
        // Define no. of descriptors
        this.DescriptorsSize = 12;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nCIC";
        this.DescriptorsNames[1] = "GATS2s";
        this.DescriptorsNames[2] = "P_VSA_m_5";
        this.DescriptorsNames[3] = "SM08_EA(bo)";
        this.DescriptorsNames[4] = "nRNNOx";
        this.DescriptorsNames[5] = "nFuranes";
        this.DescriptorsNames[6] = "N-067";
        this.DescriptorsNames[7] = "CATS2D_03_DA";
        this.DescriptorsNames[8] = "CATS2D_04_AL";
        this.DescriptorsNames[9] = "B05[O-O]";
        this.DescriptorsNames[10] = "B08[Cl-Cl]";
        this.DescriptorsNames[11] = "F04[O-Cl]";
        
        // Defines results
        this.ResultsSize = 3;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted Oral Carcinogenicity SF (log form) [log(1/(mg/kg-day))]";
        this.ResultsName[1] = "Predicted Oral Carcinogenicity SF [1/(mg/kg-day)]";
        this.ResultsName[2] = "Experimental value [1/(mg/kg-day)]";
        
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
            
            //// Custom descriptors calculation
            
            int nCIC = 0;
            int B05_O_O = 0;
            int B08_Cl_Cl = 0;
            int F04_O_Cl = 0;
            int N67 = 0;
            
            // nCIC
            nCIC = CurMolecule.GetSSSR().getAtomContainerCount();
            
            // Frequency and binary atom pairs
            int nSK = CurMolecule.GetStructure().getAtomCount();
            int[][] TopoMat = CurMolecule.GetMatrixTopologicalDistance();

            for (int i=0; i<(nSK-1); i++) {
                for (int j=i+1; j<nSK; j++) {
                    String a = CurMolecule.GetStructure().getAtom(i).getSymbol();
                    String b = CurMolecule.GetStructure().getAtom(j).getSymbol();

                    if (TopoMat[i][j] == 4) {
                        if ( ((a.equalsIgnoreCase("O")) && (b.equalsIgnoreCase("Cl"))) || ((b.equalsIgnoreCase("O")) && (a.equalsIgnoreCase("Cl"))) )
                            F04_O_Cl++;
                    }
                    if (TopoMat[i][j] == 5) {
                        if ( ((a.equalsIgnoreCase("O")) && (b.equalsIgnoreCase("O"))) )
                            B05_O_O = 1;
                    }
                    if (TopoMat[i][j] == 8) {
                        if ( ((a.equalsIgnoreCase("Cl")) && (b.equalsIgnoreCase("Cl"))) )
                            B08_Cl_Cl = 1;
                    }
                }
            }
            
            // N-67 corrected
            Pattern q1 = SmartsPattern.create("[N;D2;H;!a;!$(NC=[N,P,S,O])]([!a])[!a]", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
//            QueryAtomContainer q1 = SMARTSParser.parse("[N;D2;H;!a;!$(NC=[N,P,S,O])]([!a])[!a]");
//            CustomQueryMatcher Matcher = new CustomQueryMatcher(CurMolecule);
            boolean status = q1.matches(CurMolecule.GetStructure());
            int nmatch = 0;
//            List mappings = null;
            if (status)
                nmatch = q1.matchAll(CurMolecule.GetStructure()).countUnique();
//                mappings = Matcher.getUniqueMatchingAtoms();
//                nmatch = mappings.size();
            N67 = nmatch;
            
            
            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            Descriptors[0] = nCIC;
            Descriptors[1] = embeddedDescriptors.GATS2s;
            Descriptors[2] = embeddedDescriptors.P_VSA_m_5;
            Descriptors[3] = embeddedDescriptors.ESpm8bo;
            Descriptors[4] = embeddedDescriptors.nRNNOx;
            Descriptors[5] = embeddedDescriptors.nFuranes;
            Descriptors[6] = N67;
            Descriptors[7] = embeddedDescriptors.CATS2D_3_DA;
            Descriptors[8] = embeddedDescriptors.CATS2D_4_AL;
            Descriptors[9] = B05_O_O;
            Descriptors[10] = B08_Cl_Cl;
            Descriptors[11] = F04_O_Cl;
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }

    
    @Override
    protected short CalculateModel() {

        // Set descriptors and run ANN
        Map<String, Object> argumentsObject = new LinkedHashMap<String, Object>();
        argumentsObject.put("nCIC", Descriptors[0]);
        argumentsObject.put("GATS2s", Descriptors[1]);
        argumentsObject.put("P_VSA_m_5", Descriptors[2]);
        argumentsObject.put("SM08_EA(bo)", Descriptors[3]);
        argumentsObject.put("nRNNOx", Descriptors[4]);
        argumentsObject.put("nFuranes", Descriptors[5]);
        argumentsObject.put("N-067", Descriptors[6]);
        argumentsObject.put("CATS2D_03_DA", Descriptors[7]);
        argumentsObject.put("CATS2D_04_AL", Descriptors[8]);
        argumentsObject.put("B05[O-O]", Descriptors[9]);
        argumentsObject.put("B08[Cl-Cl]", Descriptors[10]);
        argumentsObject.put("F04[O-Cl]", Descriptors[11]);
        
        double PredictionNotNorm;
        try {
            PredictionNotNorm = ANN.Evaluate(argumentsObject);
        } catch (Exception ex) {
            return MODEL_ERROR;
        }

        // Result is normalized with the params used in the original ANN
        double Prediction = (1.0 / (0.5334846765-0.4199772985)) * PredictionNotNorm -
                (0.4199772985 / (0.5334846765-0.4199772985));
        
        CurOutput.setMainResultValue(Prediction);

        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(Prediction)); // log(1/(mg/kg-day))
        double ConvertedValue = Math.pow(10, Prediction);
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // 1/(mg/kg-day)
        else
            Res[1] = Format_4D.format(ConvertedValue); // 1/(mg/kg-day)
        Res[2] = "-";
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
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(1.5, 0.8);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(1.5, 0.8);
            ((ADIndexMaxError)CurOutput.getADIndex(ADIndexMaxError.class)).SetThresholds(1.5, 0.8);
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
        
        // Add transformed (mg/L) experimental if needed
        if (CurOutput.HasExperimental()) {
            double ConvertedValue = Math.pow(10, CurOutput.getExperimental());
            if (ConvertedValue>1)
                CurOutput.getResults()[2] = Format_2D.format(ConvertedValue); // 1/(mg/kg-day)
            else
                CurOutput.getResults()[2] = Format_4D.format(ConvertedValue); // 1/(mg/kg-day)
        }
                
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
