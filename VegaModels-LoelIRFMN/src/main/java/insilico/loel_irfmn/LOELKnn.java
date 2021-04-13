package insilico.loel_irfmn;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.GenericFailureException;
import insilico.core.knn.insilicoKnnPrediction;
import insilico.core.model.trainingset.iTrainingSet;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.similarity.SimilarMolecule;
import insilico.core.similarity.Similarity;
import insilico.core.similarity.SimilarityDescriptorsBuilder;
import lombok.extern.slf4j.Slf4j;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
@Slf4j
public class LOELKnn {
    
    private final static int NeighboursNumber = 2;
    private final static double SimThreshold = 0.85;
    private final static double SimThresholdSingle = 0.90;
    private final static double Delta = 1.0;
    
    private final iTrainingSet ModelTrainSet;
    private final SimilarityDescriptorsBuilder SimDescEngine;
    
    public LOELKnn(iTrainingSet TrainSet) {
        ModelTrainSet = TrainSet;
        SimDescEngine = new SimilarityDescriptorsBuilder();
    }

    
    public insilicoKnnPrediction Calculate(InsilicoMolecule Mol) throws GenericFailureException {
        
        // Calculate plain KNN
        insilicoKnnPrediction Res = CalculateSimpleKNN(Mol, null);
        
        // If no prediction has been provided (or an experimental value is
        // returned), return the result
        if (Res.getStatus()!= insilicoKnnPrediction.KNN_NORMAL_PREDICTION)
            return Res;
        
        // If a KNN prediction has been done, check if the neighbor molecules
        // were suitable, otherwise returns a missing value
        SimilarMolecule FirstSimMol = Res.getNeighbours().get(0);
        ArrayList<Integer> SkipMol = new ArrayList<>();
        SkipMol.add((int) FirstSimMol.getIndex());
        insilicoKnnPrediction ValidationRes = CalculateSimpleKNN(
                SmilesMolecule.Convert(ModelTrainSet.getSMILES((int) FirstSimMol.getIndex())), SkipMol);

        if (ValidationRes.getStatus() == insilicoKnnPrediction.KNN_MISSING_RANGE){
            Res.setPrediction(Descriptor.MISSING_VALUE);
            Res.setStatus(insilicoKnnPrediction.KNN_MISSING_RANGE);
            return Res;
        }
        
        if (ValidationRes.getStatus() == insilicoKnnPrediction.KNN_NORMAL_PREDICTION){
            // check if error in prediction is > 1.0
            if (Math.abs(ValidationRes.getPrediction() - ModelTrainSet.getExperimentalValue((int) FirstSimMol.getIndex())) > 1.0) {
                Res.setPrediction(Descriptor.MISSING_VALUE);
                Res.setStatus(insilicoKnnPrediction.KNN_MISSING_RANGE);
                return Res;
            }
        }
        
        return Res;
    }    

    
    public double[] CalculateOnTrainingLOO() throws GenericFailureException {
        
        double[] Pred = new double[ModelTrainSet.getMoleculesSize()];
        
        for (int i=0; i<ModelTrainSet.getMoleculesSize(); i++) {
            
            InsilicoMolecule Mol = SmilesMolecule.Convert(ModelTrainSet.getSMILES(i));
            ArrayList<Integer> SkipMol = new ArrayList<>();
            SkipMol.add(i);
            
            insilicoKnnPrediction Res = CalculateSimpleKNN(Mol, SkipMol);
            
            if (Res.getStatus()== insilicoKnnPrediction.KNN_NORMAL_PREDICTION) {

                SimilarMolecule FirstSimMol = Res.getNeighbours().get(0);
                SkipMol.add((int) FirstSimMol.getIndex());
                insilicoKnnPrediction ValidationRes = CalculateSimpleKNN(
                        SmilesMolecule.Convert(ModelTrainSet.getSMILES((int) FirstSimMol.getIndex())), SkipMol);

                if (ValidationRes.getStatus() == insilicoKnnPrediction.KNN_MISSING_RANGE){
                    Res.setPrediction(Descriptor.MISSING_VALUE);
                    Res.setStatus(insilicoKnnPrediction.KNN_MISSING_RANGE);
                }

                if (ValidationRes.getStatus() == insilicoKnnPrediction.KNN_NORMAL_PREDICTION){
                    if (Math.abs(ValidationRes.getPrediction() - ModelTrainSet.getExperimentalValue((int) FirstSimMol.getIndex())) > 1.0) {
                        Res.setPrediction(Descriptor.MISSING_VALUE);
                        Res.setStatus(insilicoKnnPrediction.KNN_MISSING_RANGE);
                    }
                }            
            }
            
            Pred[i] = Res.getPrediction();
        }
        
        return Pred;
    }    
    
    
    private insilicoKnnPrediction CalculateSimpleKNN(InsilicoMolecule Mol, ArrayList<Integer> SkipMolIndex)
            throws GenericFailureException {
        
        // Gets similar mols
        SimilarMolecule[] simmols = CalculateSimilarity(Mol, SkipMolIndex);
        
        // Get the k molecules
        ArrayList<SimilarMolecule> KNeighbours = new ArrayList<>();
        for (int i=0; i<NeighboursNumber; i++) {
            if (i >= simmols.length) 
                break;
            if (simmols[i].getSimilarity() < SimThreshold)
                break;
            
            KNeighbours.add(simmols[i]);
        }
        
        // Check if only one neighbor is available
        if (KNeighbours.size() == 1) {
            if (KNeighbours.get(0).getSimilarity() < SimThresholdSingle)
                KNeighbours = new ArrayList<>();
        }
        
        //// Possible outcomes
        
        insilicoKnnPrediction Prediction = new insilicoKnnPrediction();
        
        // 1 - No molecules
        if (KNeighbours.isEmpty()) {
            Prediction.setStatus(insilicoKnnPrediction.KNN_MISSING_NO_MOLECULES);
            return Prediction;
        }
        
        
        // 2 - Experimental found
        if (KNeighbours.get(0).getSimilarity() == 1) {
            Prediction.setPrediction(ModelTrainSet.getExperimentalValue((int)KNeighbours.get(0).getIndex()));
            Prediction.getNeighbours().add(KNeighbours.get(0));
            Prediction.setStatus(insilicoKnnPrediction.KNN_EXPERIMENTAL);
            return Prediction;
        }
        
        
        // 3 - Check experimental delta (if k>1)
        if (KNeighbours.size() > 1) {
            double min = ModelTrainSet.getExperimentalValue((int)KNeighbours.get(0).getIndex());
            double max = min;
            for (SimilarMolecule curMol : KNeighbours) {
                double curNeighExp = ModelTrainSet.getExperimentalValue((int)curMol.getIndex());
                min = Math.min(curNeighExp, min);
                max = Math.max(curNeighExp, max);
            }

            if ( Math.abs(max-min) > Delta ) {
                Prediction.setStatus(insilicoKnnPrediction.KNN_MISSING_RANGE);
                Prediction.setNeighbours(KNeighbours);
                return Prediction;
            }
        }

        // 4 - Prediction on one or more molecules
        double Pred = 0.0;
        for (SimilarMolecule curMol : KNeighbours)
            Pred += ModelTrainSet.getExperimentalValue((int)curMol.getIndex());
        Pred /= (double)KNeighbours.size();

        Prediction.setPrediction(Pred);
        Prediction.setNeighbours(KNeighbours);
        Prediction.setStatus(insilicoKnnPrediction.KNN_NORMAL_PREDICTION);
        return Prediction;
    }
    
    
    private SimilarMolecule[] CalculateSimilarity(InsilicoMolecule Mol,
                                                  ArrayList<Integer> SkipMolIndex) throws GenericFailureException {
        
        SimilarMolecule[] SimilarMols;
        
        if ((ModelTrainSet == null)||(ModelTrainSet.getMoleculesSize()==0)) {
            SimilarMols = null;
            log.warn("Unable to retrieve training set for similarity calculation");
            throw new GenericFailureException("Unable to retrieve training set");
        }
        
        Similarity SIM = new Similarity();
        SimilarMols = new SimilarMolecule[ModelTrainSet.getMoleculesSize()];
        
        // Calculates similarity for all molecules
        for (int idx=0; idx<ModelTrainSet.getMoleculesSize(); idx++) {
            
            // just for LOO execution, skip calculation for the molecule
            // with index SkipMolIndex (set to null for normal execution)
            boolean SkipMolecule = false;
            if (SkipMolIndex != null)
                for (Integer skipIdx : SkipMolIndex)
                    if (idx == skipIdx) {
                        SkipMolecule = true;
                        break;
                    }
            if (SkipMolecule) {
                SimilarMols[idx] = new SimilarMolecule(idx, 0.0);
                continue;
            }
            
            if (!Mol.HasSimilarityDescriptors()){

                Mol.SetSimilarityDescriptors(SimDescEngine.Calculate(Mol));
            }
            double curSim;
            try {
                curSim = SIM.Calculate(Mol.GetSimilarityDescriptors(),
                        ModelTrainSet.getSimilarityDescriptor(idx));
                
                // If similarity is equal to 1, checks if mols are really identical
                if (curSim == 1.0) {
                    IAtomContainer A = Mol.GetStructure();
                    IAtomContainer B = SmilesMolecule.Convert(ModelTrainSet.getSMILES(idx)).GetStructure();
                    boolean AreIsomorph = Similarity.CheckIsomorphism(A, B);
                    if (!AreIsomorph)
                        curSim = 0.999;
                }

                // Adjust low similarities to avoid problems in some following indices
                if (curSim < 0.38)
                    curSim = 0.38;
                
            } catch (Throwable e) {
                log.warn("Similarity calculation: unable to calculate for training set molecule "
                        + idx + ": " + ModelTrainSet.getSMILES(idx));
                curSim = 0;
            }
            SimilarMols[idx] = new SimilarMolecule(idx, curSim);
        }

        // Sorts the array of similar molecules (index=0 means most similar molecule)
        Arrays.sort(SimilarMols);
        return SimilarMols;
    }     
    
    
}
