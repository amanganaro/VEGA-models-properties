package insilico.fish_nic;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQuantitative;
import insilico.core.ad.item.*;
import insilico.core.constant.MessagesAD;
import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.descriptor.blocks.*;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;
import insilico.core.molecule.tools.Manipulator;
import insilico.core.tools.utils.ModelUtilities;
import insilico.fish_nic.descriptors.EmbeddedDescriptors;
import insilico.fish_nic.descriptors.weights.CovalentRadius;
import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class ismFishNic extends InsilicoModel {
    
    private static final long serialVersionUID = 1L;
    
    private static final String ModelData = "/data/model_fish_nic.xml";
    
    private final modelCPANN CPANN;
//    private final ADNeuronsChecker NeuronChecker;
    private double[] ANNResults;
    private double MW;
    
    
    public ismFishNic() 
            throws InitFailureException {
        super(ModelData);
        
        // Builds model objects
        CPANN = new modelCPANN();
//        NeuronChecker = new ADNeuronsChecker(35, 35, "/insilico/model/carcinogenicity/caesar/data/neurons.txt");

        // Define no. of descriptors
        this.DescriptorsSize = 20;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "nCIR";
        this.DescriptorsNames[1] = "nTB";
        this.DescriptorsNames[2] = "nP";
        this.DescriptorsNames[3] = "nCL";
        this.DescriptorsNames[4] = "nR10";
        this.DescriptorsNames[5] = "TI1";
        this.DescriptorsNames[6] = "S2K";
        this.DescriptorsNames[7] = "T(N..P)";
        this.DescriptorsNames[8] = "T(P..Cl)";
        this.DescriptorsNames[9] = "T(Cl..Cl)";
        this.DescriptorsNames[10] = "piPC09";
        this.DescriptorsNames[11] = "PCR";
        this.DescriptorsNames[12] = "Xindex";
        this.DescriptorsNames[13] = "MATS1e";
        this.DescriptorsNames[14] = "GATS7m";
        this.DescriptorsNames[15] = "EEig14x";
        this.DescriptorsNames[16] = "EEig14d";
        this.DescriptorsNames[17] = "ESpm01d";
        this.DescriptorsNames[18] = "GGI4";
        this.DescriptorsNames[19] = "Seigv";

        
        // Defines results
        this.ResultsSize = 4;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted fish toxicity [log(1/(mmol/L))]";
        this.ResultsName[1] = "Predicted fish toxicity [mg/l]";
        this.ResultsName[2] = "Molecular Weight";
        this.ResultsName[3] = "Experimental value [mg/l]";
        
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
            
            // Custom descriptors
            // TO DO : move in proper descriptors blocks

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);
            
            int nSK = CurMolecule.GetStructure().getAtomCount();
            int nBO = CurMolecule.GetStructure().getBondCount();

            
            ////// nCIR
            
            int nCIR = CurMolecule.GetAllRings().getAtomContainerCount();
            
            
            ////// TI1
            
            int[][] mat = CurMolecule.GetMatrixLaplace();
            
            double[][] eig_mat = new double[mat.length][mat[0].length];
            for (int k=0; k<mat.length; k++)
                for (int z=0; z<mat[0].length; z++)
                    eig_mat[k][z] = mat[k][z];
            
            // Calculates eigenvalues
            Matrix DataMatrix = new Matrix(eig_mat);
            double[] eigenvalues;
            EigenvalueDecomposition ed = new EigenvalueDecomposition(DataMatrix);
            eigenvalues = ed.getRealEigenvalues();
            //Arrays.sort(eigenvalues);
            
            double eig_sum = 0;
            for (double eig : eigenvalues)
                if (eig>0.000001)
                    eig_sum += 1.0/eig;
            
            double QW = eig_sum * (double)nSK;
            double TI1 = 2 * QW * Math.log10((double)nBO/(double)nSK);

            
            ////// S2K
            
            double S2K = Descriptor.MISSING_VALUE;
            
            // TotPC2
            double TotPC2 = embeddedDescriptors.getMPC2();
            
            // Alpha
            double[] KierCR = getKierCovRadius(CurMolecule.GetStructure());
            boolean SomeMissing = false;
            for (double d : KierCR)
                if (d== Descriptor.MISSING_VALUE)
                    SomeMissing = true;
            
            if (!SomeMissing) {
                double Csp3_radius =  0.77;
                double alpha = 0;
                for (double d : KierCR)
                    alpha += (d / Csp3_radius) - 1;

                // S2K
                if (nSK == 1)
                    S2K = 0;
                else if (nSK == 2)
                    S2K = 1 +alpha;
                else {
                    if (TotPC2 == 0)
                        S2K = 0;
                    else
                        S2K = (nSK + alpha - 1) * Math.pow((nSK + alpha - 2),2) /
                                Math.pow((TotPC2 + alpha),2);
                }
            } else {
                S2K = Descriptor.MISSING_VALUE;
            }
            
            
            ////// Xindex
            
            int nEdges = CurMolecule.GetStructure().getBondCount();
            int nH = 0;
            for (IAtom a : CurMolecule.GetStructure().atoms()) {
                int H = Manipulator.CountImplicitHydrogens( a);
                nH += H;
                nEdges += H;
            }
            
            int[][] DistanceT = CurMolecule.GetMatrixTopologicalDistance();
            
            // Calculate distance degree and distance frequencies for all atoms
            int[] DistDeg = new int[nSK];
            int [][] AtomDistanceFreq = new int[nSK][nSK+1];

            for (int at=0; at<nSK; at++) {
                DistDeg[at] = 0;
                for (int adf : AtomDistanceFreq[at]) adf = 0;
                
                for (int z=0; z<nSK; z++) {
                    if (at==z) continue;
                    DistDeg[at] += DistanceT[at][z];
                    AtomDistanceFreq[at][DistanceT[at][z]]++;
                }
            }
            
            // Calculate Xlocal (weight for Xindex)
            double[] Xlocal = new double[nSK];
            double[] Ylocal = new double[nSK];
            
            for (int at=0; at<nSK; at++) {
                Ylocal[at] = 0;
                int idx=0;
                while ( (idx<nSK) && (AtomDistanceFreq[at][idx+1]>0) ) {
                    double val = AtomDistanceFreq[at][idx+1] * (idx+1);
                    Ylocal[at] += val * Math.log(idx+1) * (1.0 / Math.log(2));
                    idx++;
                }
            }            

            for (int at=0; at<nSK; at++) {
                double Sig = DistDeg[at] * Math.log(DistDeg[at]) * (1.0 / Math.log(2));
                Xlocal[at] = Sig - Ylocal[at];
            }

            // Calculate Xindex
            double Xindex = 0;
            for (int at=0; at<nSK; at++) 
                for (int at2=0; at2<nSK; at2++)
                    if (DistanceT[at][at2] == 1) 
                        if (Xlocal[at] > 0)
                            Xindex += 1.0 / Math.sqrt(Xlocal[at] * Xlocal[at2]);
            Xindex = Xindex / 2.0 * (nEdges - (nSK+nH) + nSK ) / (nEdges - (nSK+nH) + 2);
          
              
            Descriptors = new double[DescriptorsSize];
            Descriptors[0] = nCIR;
            Descriptors[1] = embeddedDescriptors.getNTB();
            Descriptors[2] = embeddedDescriptors.getNP();
            Descriptors[3] = embeddedDescriptors.getNCL();
            Descriptors[4] = embeddedDescriptors.getNR10();
            Descriptors[5] = TI1;
            Descriptors[6] = S2K;
            Descriptors[7] = embeddedDescriptors.getT_N_P();
            Descriptors[8] = embeddedDescriptors.getT_P_CL();
            Descriptors[9] = embeddedDescriptors.getT_CL_CL();
            Descriptors[10] = embeddedDescriptors.getPiPC09();
            Descriptors[11] = embeddedDescriptors.getPCR();
            Descriptors[12] = Xindex;
            Descriptors[13] = embeddedDescriptors.getMATS1e();
            Descriptors[14] = embeddedDescriptors.getGATS7m();
            Descriptors[15] = embeddedDescriptors.getEEig14ed(); // non va
            Descriptors[16] = embeddedDescriptors.getEEig14dm(); // non va
            Descriptors[17] = embeddedDescriptors.getESpm1dm();
            Descriptors[18] = embeddedDescriptors.getGGI4();
            Descriptors[19] = embeddedDescriptors.getSeigv();
                        
            // MW in constitutional is given as a SCALED 
            // value (on carbon). Here it is transformed in real values
            double CarbonWeight = 12.011;
            MW = CarbonWeight * embeddedDescriptors.getMw();
            
        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }
    
        return DESCRIPTORS_CALCULATED;
    }


    
    @Override
    protected short CalculateModel() {
        
        ANNResults = CPANN.CalculatePrediction(Descriptors);
        
        double MainResult = ANNResults[0];
        CurOutput.setMainResultValue(MainResult);
        
        String[] Res = new String[ResultsSize];
        Res[0] = String.valueOf(Format_2D.format(MainResult)); // toxicity log(1/(mol/L))
        double ConvertedValue = Math.pow(10, (-1 * MainResult)) * MW;
        if (ConvertedValue>1)
            Res[1] = Format_2D.format(ConvertedValue); // in mg/L
        else
            Res[1] = Format_4D.format(ConvertedValue); // in mg/L
        Res[2] = Format_2D.format(MW); // MW
        Res[3] = "-"; // Converted experimental - set after
        
        CurOutput.setResults(Res);
        
        return MODEL_CALCULATED;
    }

    
    @Override
    protected short CalculateAD() {

        //// TODO indici relativi al nn
        
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
        
        // Add transformed (mg/L) experimental if needed
        if (CurOutput.HasExperimental()) {
            double ConvertedExp = Math.pow(10, (-1 * CurOutput.getExperimental())) * MW;
            if (ConvertedExp>1)
                CurOutput.getResults()[3] = Format_2D.format(ConvertedExp); // in mg/L
            else
                CurOutput.getResults()[3] = Format_4D.format(ConvertedExp); // in mg/L
        }
                        
        return InsilicoModel.AD_CALCULATED;
    }
        

    @Override
    protected void CalculateAssessment() {
        
        double LC_threshold_red = 1; // in mg/l
        double LC_threshold_orange = 10; // in mg/l
        double LC_threshold_yellow = 100; // in mg/l
        
        // Sets assessment message
        // Can't use default utilities because a different experimental has
        // to be set (mg/L) if available

        String ADItemWarnings =
                ModelUtilities.BuildADItemsWarningMsg(CurOutput.getADIndex());
        
        String Result = CurOutput.getResults()[1] + " mg/L";
        
        switch (CurOutput.getADI().GetAssessmentClass()) {
            case ADIndex.INDEX_LOW:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_LOW, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_LOW, Result, ADItemWarnings));
                break;
            case ADIndex.INDEX_MEDIUM:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_MEDIUM, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_MEDIUM, Result, ADItemWarnings));
                break;
            case ADIndex.INDEX_HIGH:
                CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_HIGH, Result));
                CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_HIGH, Result));
                if (!ADItemWarnings.isEmpty())
                    CurOutput.setAssessmentVerbose(CurOutput.getAssessmentVerbose() +
                            String.format(MessagesAD.ASSESS_LONG_ADD_ISSUES, ADItemWarnings));
                break;
        }

        // Override assessment if experimental value is available
        if (CurOutput.HasExperimental()) {
            CurOutput.setAssessmentVerbose(String.format(MessagesAD.ASSESS_LONG_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L", CurOutput.getAssessment()));
            CurOutput.setAssessment(String.format(MessagesAD.ASSESS_SHORT_EXPERIMENTAL, CurOutput.getResults()[3] + " mg/L"));
        }
        

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        Val = Math.pow(10, (-1 * Val)) * MW; // convert to mg/L
        if (Val < LC_threshold_red)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else if (Val < LC_threshold_orange)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_ORANGE);
        else if (Val < LC_threshold_yellow)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        
    }
    
    
    private static double[] getKierCovRadius(IAtomContainer mol) {
        
        int nSK = mol.getAtomCount();
        double[] covrad = CovalentRadius.getWeights(mol);
        
        // Calculate Kier radius for specific atom types, if not found normal
        // covalent radius is used
        
        for (int i=0; i<nSK; i++) {
            IAtom at = mol.getAtom(i);
            double rad = Descriptor.MISSING_VALUE;
//            int vd = mol.getConnectedAtomsCount(at) + Manipulator.CountImplicitHydrogens(at);
            int vd = mol.getConnectedBondsCount(at);
            
            if (at.getSymbol().compareTo("C")==0) {
                if (vd==4) rad = 0.77; // C sp3 
                if (vd==3) rad = 0.67; // C sp2 
                if (vd==2) rad = 0.60; // C sp 
                if (vd==1) rad = 0.0; // 
            }
            if (at.getSymbol().compareTo("N")==0) {
                if (vd>2) rad = 0.74; // N sp3
                if (vd==2) rad = 0.62; // N sp2
                if (vd==1) rad = 0.55; // N sp
            }
            if (at.getSymbol().compareTo("O")==0) {
                if (vd>1) rad = 0.74; // O sp3
                if (vd==1) rad = 0.62; // O sp2
            }
            if (at.getSymbol().compareTo("S")==0) {
                if (vd>1) rad = 1.04; // O sp3
                if (vd==1) rad = 0.94; // O sp2
            }
            if (at.getSymbol().compareTo("P")==0) {
                if (vd>2) rad = 1.10; // P sp3
                if (vd==2) rad = 1.00; // P sp2
            }
            if (at.getSymbol().compareTo("F")==0)
                rad = 0.72;
            if (at.getSymbol().compareTo("Cl")==0)
                rad = 0.99;
            if (at.getSymbol().compareTo("Br")==0)
                rad = 1.14;
            if (at.getSymbol().compareTo("I")==0)
                rad = 1.33;
            
            if (rad != Descriptor.MISSING_VALUE)
                covrad[i] = rad;
        }
        
        return covrad;
    }    
}
