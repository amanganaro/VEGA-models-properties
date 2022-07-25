package insilico.devtox_caesar;

import insilico.core.ad.ADCheckACF;
import insilico.core.ad.ADCheckDescriptorRange;
import insilico.core.ad.ADCheckIndicesQualitative;
import insilico.core.ad.item.*;
import insilico.core.descriptor.DescriptorsEngine;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.core.model.InsilicoModelOutput;

import insilico.core.tools.utils.ModelUtilities;
import insilico.devtox_caesar.descriptors.EmbeddedDescriptors;
import insilico.devtox_caesar.descriptors.IcycemDescriptor;
import lombok.extern.log4j.Log4j;

import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
@Log4j
public class ismDevtoxCaesar extends InsilicoModel {

    private static final long serialVersionUID = 1L;

    private static final String ModelData = "/data/model_devtox_caesar.xml";

    private final modelRandomForest RandomForest;


    public ismDevtoxCaesar()
            throws InitFailureException {
        super(ModelData);

        // Builds model objects
        RandomForest = new modelRandomForest("/data/sel18_EPA_training.arff");

        // Define no. of descriptors
        this.DescriptorsSize = 13;
        this.DescriptorsNames = new String[DescriptorsSize];
        this.DescriptorsNames[0] = "icycem";
        this.DescriptorsNames[1] = "BEHm1";
        this.DescriptorsNames[2] = "BELp3";
        this.DescriptorsNames[3] = "BELv1";
        this.DescriptorsNames[4] = "BELv8";
        this.DescriptorsNames[5] = "GATS1p";
        this.DescriptorsNames[6] = "GATS2m";
        this.DescriptorsNames[7] = "GATS3v";
        this.DescriptorsNames[8] = "MATS1p";
        this.DescriptorsNames[9] = "MATS4p";
        this.DescriptorsNames[10] = "MATS4v";
        this.DescriptorsNames[11] = "SdssC";
        this.DescriptorsNames[12] = "SHssNH";

        // Defines results
        this.ResultsSize = 1;
        this.ResultsName = new String[ResultsSize];
        this.ResultsName[0] = "Predicted developmental toxicity activity";

        // Define AD items
        this.ADItemsName = new String[5];
        this.ADItemsName[0] = new ADIndexSimilarity().GetIndexName();
        this.ADItemsName[1] = new ADIndexAccuracy().GetIndexName();
        this.ADItemsName[2] = new ADIndexConcordance().GetIndexName();
        this.ADItemsName[3] = new ADIndexRange().GetIndexName();
        this.ADItemsName[4] = new ADIndexACF().GetIndexName();

    }


    @Override
    protected short CalculateDescriptors(DescriptorsEngine DescEngine) {

        try {

            Descriptors = new double[DescriptorsSize];

            EmbeddedDescriptors embeddedDescriptors = new EmbeddedDescriptors(CurMolecule);

            // icycem is calculated with a custom class
            Descriptors[0] = IcycemDescriptor.Calculate(CurMolecule);

            Descriptors[1] = embeddedDescriptors.getBEH1m();

            // Burden low eigenvalue descriptors (BEL) are different from
            // the descriptors calculated in the original custom code,
            // they just have opposite sign
            Descriptors[2] = (-1) * embeddedDescriptors.getBEL3p();
            Descriptors[3] = (-1) * embeddedDescriptors.getBEL1v();
            Descriptors[4] = (-1) * embeddedDescriptors.getBEL8v();

            Descriptors[5] = embeddedDescriptors.getGATS1p();
            Descriptors[6] = embeddedDescriptors.getGATS2m();
            Descriptors[7] = embeddedDescriptors.getGATS3v();
            Descriptors[8] = embeddedDescriptors.getMATS1p();
            Descriptors[9] = embeddedDescriptors.getMATS4p();
            Descriptors[10] = embeddedDescriptors.getMATS4v();
            Descriptors[11] = embeddedDescriptors.getSdssC();
            Descriptors[12] = embeddedDescriptors.getSHssNH();

        } catch (Throwable e) {
            return DESCRIPTORS_ERROR;
        }

        return DESCRIPTORS_CALCULATED;
    }



    @Override
    protected short CalculateModel() {

        int MainResult = -1;

        try {

            // Creates the ARFF file
            String header = "@relation model_test"+"\n"+"\n"+"@attribute SdssC numeric"+
                "\n"+"@attribute SHssNH numeric"+"\n"+"@attribute icycem" +
                " numeric"+"\n"+"@attribute BEHm1 numeric"+"\n"+"@attribute " +
                "BELv1 numeric"+"\n"+"@attribute BELv8 numeric"+"\n"+
                "@attribute" +" BELp3 numeric"+"\n"+"@attribute MATS4v numeric"
                +"\n"+"@attribute MATS1p numeric"+"\n"+"@attribute MATS4p" +
                " numeric"+"\n"+"@attribute GATS2m numeric"+"\n"+"@attribute" +
                " GATS3v" +" numeric"+"\n"+"@attribute GATS1p numeric"+"\n"+
                "@attribute" +" Tox {N,D}"+"\n"+"\n"+"@data"+"\n";

            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            OutputStreamWriter osw = new OutputStreamWriter(baos);
            BufferedWriter bw = new BufferedWriter(osw);

            bw.write(header);
            bw.write(Double.toString(Descriptors[11]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[12]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[0]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[1]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[3]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[4]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[2]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[10]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[8]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[9]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[6]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[7]));
            bw.write(",");
            bw.write(Double.toString(Descriptors[5]));
            bw.write(",");
            bw.write("?");

            bw.close();

            byte[] Input_ARFF = baos.toByteArray();

            boolean CurResult = RandomForest.ExecuteModel(Input_ARFF);
            if (CurResult)
                MainResult = 1;
            else
                MainResult = 0;

            // No undefined class can be returned

        } catch (IOException | GenericFailureException e) {
            return MODEL_ERROR;
        }

        CurOutput.setMainResultValue(MainResult);

        String[] Res = new String[ResultsSize];
        try {
            Res[0] = this.GetTrainingSet().getClassLabel(MainResult); // devtox classification
        } catch (Throwable ex) {
            log.warn("Unable to find label for toxicity value " + MainResult);
            Res[0] = Integer.toString(MainResult);
        }

        CurOutput.setResults(Res);

        return MODEL_CALCULATED;
    }


    @Override
    protected short CalculateAD() {

        // Calculates various AD indices
        ADCheckIndicesQualitative adq = new ADCheckIndicesQualitative(TS);
        adq.AddMappingToPositiveValue(1);
        adq.AddMappingToNegativeValue(0);
        adq.setMoleculesForIndexSize(2);
        if (!adq.Calculate(CurMolecule, CurOutput))
            return InsilicoModel.AD_ERROR;

        // Sets threshold for AD indices
        try {
            ((ADIndexSimilarity)CurOutput.getADIndex(ADIndexSimilarity.class)).SetThresholds(0.8, 0.7);
            ((ADIndexAccuracy)CurOutput.getADIndex(ADIndexAccuracy.class)).SetThresholds(0.8, 0.6);
            ((ADIndexConcordance)CurOutput.getADIndex(ADIndexConcordance.class)).SetThresholds(0.8, 0.6);
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

        ADIndexADI ADI = new ADIndexADI();
        ADI.SetIndexValue(ADIValue);
        ADI.SetThresholds(0.8, 0.7);
        CurOutput.setADI(ADI);

        return InsilicoModel.AD_CALCULATED;
    }


    @Override
    protected void CalculateAssessment() {

        // Sets assessment message
        ModelUtilities.SetDefaultAssessment(CurOutput, CurOutput.getResults()[0]);

        // Sets assessment status
        double Val = CurOutput.HasExperimental() ? CurOutput.getExperimental() : CurOutput.getMainResultValue();
        if (Val == 0)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_GREEN);
        else if (Val == 1)
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_RED);
        else
            CurOutput.setAssessmentStatus(InsilicoModelOutput.ASSESS_YELLOW);

    }

}
