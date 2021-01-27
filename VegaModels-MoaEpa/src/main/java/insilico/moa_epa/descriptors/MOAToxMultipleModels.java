package insilico.moa_epa.descriptors;

/**
 * 
 * Implementation for the LDA models for the 7 MOA classes, as described
 * in the original publication used as reference for the T.E.S.T. software
 * 
 * MANCANO:
 * (IC) ok non torna qualche val
 * xch10 (CHI chain 10)
 * 
 * 
 * 
 * @author Alberto
 */
public class MOAToxMultipleModels {

    // coeff for single MOA models
    
    private final static double[] COEFF_AChE_Inhibition = {
        + 0.8416, // (SdsssP_acnt)
        - 0.0055, // (ic)
        - 0.0033, // (MDEC22)
        - 0.0779, // (MATS1m)
        + 0.2968, // (MATS2e)
        + 0.2071, // (GATS2e)
        - 0.7224, // (-O- [2 phosphorus attach])
        - 0.4428, // (-OH [phosphorus attach])
        + 0.0745, // (-NH- [aliphatic attach])
        - 0.1102, // (-NH- [aromatic attach])
        + 0.5166, // (-C(=O)O- [nitrogen attach])
        - 0.1102    
    };
    
    private final static double[] COEFF_Anticoagulation = {
        + 13.1706, // (xch10)
        - 0.0098, // (MDEC13)
        + 0.0025, // (MDEC23)
        - 0.0788, // (BELe1)
        - 0.0963, // (nR10)
        + 0.0164, // (GATS6v)
        + 0.0870, // (-CH< [aromatic attach])
        + 0.1001, // (>C= [aromatic attach])
        + 0.1896, // (C=O(ketone, aliphatic attach))
        + 0.1836, // (-C(=O)- [aromatic attach])
        + 0.2191, // (-C(=O)O- [cyclic])
        + 0.1148    
    };
    
    private final static double[] COEFF_nACHr_Antagonism = {
        + 0.0020, // (Gmin)
        + 0.6646, // (-C#N [aliphatic nitrogen attach])
        + 0.9979, // (-NO2 [nitrogen attach])
        + 0.0035        
    };
    
    private final static double[] COEFF_Narcosis = {
        + 0.3290, // (SdssNp)
        - 0.3621, // (SdCH2_acnt)
        - 0.2639, // (StsC_acnt)
        - 0.7252, // (SdsssP_acnt)
        - 0.0231, // (MDEC33)
        - 0.0157, // (AMW)
        - 0.0007, // (SRW07)
        - 0.1787, // (-CH< [aromatic attach])
        - 0.7161, // (-NH- [nitrogen attach])
        - 0.9070, // (-CHO [aliphatic attach])
        - 0.6308, // (-C(=O)O- [nitrogen attach])
        + 1.0007        
    };
    
    private final static double[] COEFF_Neurotoxicity = {
        + 0.0614, // (SdssC)
        + 0.0115, // (StN)
        + 0.0054, // (ic)
        - 0.0935, // (icycem)
        + 0.0654, // (MDEC34)
        + 0.0103, // (nX)
        - 0.0628, // (nR05)
        + 0.0958, // (nR06)
        - 0.1398, // (nR10)
        + 0.0008, // (SRW07)
        + 0.1499, // (-C(=O)O- [aliphatic attach])
        - 0.0477        
    };
    
    private final static double[] COEFF_Reactivity = {
        + 0.0679, // (StsC)
        + 0.0120, // (SsOm)
        + 0.3471, // (SdCH2_acnt)
        + 0.2950, // (StCH_acnt)
        - 0.1039, // (BEHp4)
        + 0.5010, // (-NH- [nitrogen attach])
        + 0.3237, // (-S- [sulfur attach])
        + 0.1813, // (-C(=O)- [aromatic attach])
        + 0.8807, // (-CHO [aliphatic attach])
        + 0.6064, // (-CHO [aromatic attach])
        + 0.5932, // (CH2=CHC(=O)O-)
        + 0.2951        
    };
    
    private final static double[] COEFF_Uncoupler = {
        - 0.1872, // (SdssNp)
        + 0.0751, // (Hmin)
        + 0.2243, // (Qsv)
        - 0.0962, // (ivdem)
        + 0.1385, // (Ms)
        + 0.0388, // (ATS4m)
        - 0.0271, // (GATS8m)
        + 0.1425, // (-Br [aromatic attach])
        + 0.0389, // (-Cl [aromatic attach])
        + 0.3286, // (-I [aromatic attach])
        + 0.2637, // (-OH [aromatic attach])
        - 0.5566        
    };
    
    
    public final static int DescriptorsSize = 11 + 11 + 3 + 11 + 11 + 11 + 11;

    public final static int MOAModelsNumber = 7;
    public final static String[] MOANames = {
        "AChE Inhibition",
        "Anticoagulation",
        "nACHr Antagonism",
        "Narcosis",
        "Neurotoxicity",
        "Reactivity",
        "Uncoupler"
    };
    
        
    /**
     * 
     * @param Descriptors
     * @return 
     */
    public static double[] CalculateScores(double[] Descriptors) {
        
        //      Input descriptors: 
        //
        // 	SdsssP_acnt
        // 	ic
        // 	MDEC22
        // 	MATS1m
        // 	MATS2e
        // 	GATS2e
        // 	-O- [2 phosphorus attach]
        // 	-OH [phosphorus attach]
        // 	-NH- [aliphatic attach]
        // 	-NH- [aromatic attach]
        // 	-C(=O)O- [nitrogen attach]
        // 	xch10
        // 	MDEC13
        // 	MDEC23
        // 	BELe1
        // 	nR10
        // 	GATS6v
        // 	-CH< [aromatic attach]
        // 	>C= [aromatic attach]
        // 	C=O(ketone, aliphatic attach)
        // 	-C(=O)- [aromatic attach]
        // 	-C(=O)O- [cyclic]
        // 	Gmin
        // 	-C#N [aliphatic nitrogen attach]
        // 	-NO2 [nitrogen attach]
        // 	SdssNp
        // 	SdCH2_acnt
        // 	StsC_acnt
        // 	SdsssP_acnt
        // 	MDEC33
        // 	AMW
        // 	SRW07
        // 	-CH< [aromatic attach]
        // 	-NH- [nitrogen attach]
        // 	-CHO [aliphatic attach]
        // 	-C(=O)O- [nitrogen attach]
        // 	SdssC
        // 	StN
        // 	ic
        // 	icycem
        // 	MDEC34
        // 	nX
        // 	nR05
        // 	nR06
        // 	nR10
        // 	SRW07
        // 	-C(=O)O- [aliphatic attach]
        // 	StsC
        // 	SsOm
        // 	SdCH2_acnt
        // 	StCH_acnt
        // 	BEHp4
        // 	-NH- [nitrogen attach]
        // 	-S- [sulfur attach]
        // 	-C(=O)- [aromatic attach]
        // 	-CHO [aliphatic attach]
        // 	-CHO [aromatic attach]
        // 	CH2=CHC(=O)O-
        // 	SdssNp
        // 	Hmin
        // 	Qsv
        // 	ivdem
        // 	Ms
        // 	ATS4m
        // 	GATS8m
        // 	-Br [aromatic attach]
        // 	-Cl [aromatic attach]
        // 	-I [aromatic attach]
        // 	-OH [aromatic attach]        
        
        double[] MOAScores = new double[MOAModelsNumber];
        
        //// 1. AChE Inhibition
        MOAScores[0] = 0;
        for (int i=0; i<11; i++)
            MOAScores[0] += Descriptors[i] * COEFF_AChE_Inhibition[i];
        MOAScores[0] += COEFF_AChE_Inhibition[11]; // intercept


        //// 2. Anticoagulation
        MOAScores[1] = 0;
        for (int i=0; i<11; i++)
            MOAScores[1] += Descriptors[11 + i] * COEFF_Anticoagulation[i];
        MOAScores[1] += COEFF_Anticoagulation[11]; // intercept


        //// 3. nACHr Antagonism
        MOAScores[2] = 0;
        for (int i=0; i<3; i++)
            MOAScores[2] += Descriptors[22 + i] * COEFF_nACHr_Antagonism[i];
        MOAScores[2] += COEFF_nACHr_Antagonism[3]; // intercept


        //// 4. Narcosis
        MOAScores[3] = 0;
        for (int i=0; i<11; i++)
            MOAScores[3] += Descriptors[25 + i] * COEFF_Narcosis[i];
        MOAScores[3] += COEFF_Narcosis[11]; // intercept


        //// 5. Neurotoxicity
        MOAScores[4] = 0;
        for (int i=0; i<11; i++)
            MOAScores[4] += Descriptors[36 + i] * COEFF_Neurotoxicity[i];
        MOAScores[4] += COEFF_Neurotoxicity[11]; // intercept


        //// 6. Reactivity
        MOAScores[5] = 0;
        for (int i=0; i<11; i++)
            MOAScores[5] += Descriptors[47 + i] * COEFF_Reactivity[i];
        MOAScores[5] += COEFF_Reactivity[11]; // intercept


        //// 7. Uncoupler
        MOAScores[6] = 0;
        for (int i=0; i<11; i++)
            MOAScores[6] += Descriptors[58 + i] * COEFF_Uncoupler[i];
        MOAScores[6] += COEFF_Uncoupler[11]; // intercept        
        
        return MOAScores;
    }
    
    
    // test dataset taken from the T.E.S.T. software
    
    public static String[] MOADataset = {
        "O=C2NC(=O)C(c1ccccc1)(C(=O)N2)CC",
        "c1cc(ccc1C(c2ccc(cc2)Cl)C(Cl)(Cl)Cl)Cl",
        "O=C(O)c1c(ccc(c1Cl)Cl)Cl",
        "O=[N+]([O-])c1ccc(O)c(c1)[N+](=O)[O-]",
        "O=C(OCC)N",
        "O=P(OC)(OC)C(O)C(Cl)(Cl)Cl",
        "O=NN(CC)CC",
        "O=C(N)c1ccccc1",
        "O(c1ccc(c(c1)C)SC)P(OC)(OC)=S",
        "O=[N+]([O-])OCC(O[N+](=O)[O-])CO[N+](=O)[O-]",
        "C(Cl)(Cl)(Cl)Cl",
        "O=[N+]([O-])c1ccc(OP(OCC)(OCC)=S)cc1",
        "C(=NCC=C)=S",
        "NN(C)C",
        "OC(C)(C)C(Cl)(Cl)Cl",
        "O=C1NC(=O)C(C(=O)N1)(CC)CCC(C)C",
        "OCC(O)C",
        "C2(=C(C3(C1C(CC(C1Cl)Cl)C2(C3(Cl)Cl)Cl)Cl)Cl)Cl",
        "O=C1C=C(C(=O)c2ccccc12)C",
        "C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)Cl",
        "Oc1c(cc(c(c1Cl)Cl)Cl)Cl",
        "Oc1ccc(c(c1)C)Cl",
        "O=C(O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O",
        "O(CC)CC",
        "O1C2C1C3CC2C4C3C5(C(=C(C4(C5(Cl)Cl)Cl)Cl)Cl)Cl",
        "Nc1ccccc1",
        "NC(C)=S",
        "O=P(OC=C(Cl)Cl)(OC)OC",
        "O=NN(C)C",
        "O=C(Oc1cccc2ccccc12)NC",
        "OCC",
        "O=C(O)C",
        "O=C(N)c1ccccc1(O)",
        "O=CCCCCC",
        "O=C1Oc4ccccc4(C(O)=C1CC=2C(=O)Oc3ccccc3(C=2(O)))",
        "O=Cc2ccc(Oc1ccccc1)cc2",
        "OC",
        "OC(C)C",
        "O=C(C)C",
        "C(Cl)(Cl)Cl",
        "C(C(Cl)(Cl)Cl)(Cl)(Cl)Cl",
        "O=C(O)CS",
        "Oc1c(cc(c(c1Cc2c(O)c(cc(c2Cl)Cl)Cl)Cl)Cl)Cl",
        "O=C(OCc1ccc(cc1C)C)C2C(C=C(C)C)C2(C)(C)",
        "O=C(c1ccc(N)cc1)CC",
        "OCCC",
        "OCCCC",
        "OCCCCC",
        "c1ccccc1",
        "CC(Cl)(Cl)Cl",
        "O(c1ccc(cc1)C(c2ccc(OC)cc2)C(Cl)(Cl)Cl)C",
        "c1cc(ccc1C(c2ccc(cc2)Cl)C(Cl)Cl)Cl",
        "N#CC",
        "O=CC",
        "C(Cl)Cl",
        "O1CC1",
        "C=C(Cl)Cl",
        "O=C(C)Cl",
        "C(I)(I)I",
        "OC(C)(C)C",
        "FC(F)(F)CO",
        "O=C(C)C(C)(C)C",
        "O=C(O)C(C)(Cl)Cl",
        "C(C(Cl)(Cl)Cl)(Cl)Cl",
        "O=C(O)C(Cl)(Cl)Cl",
        "O=C1C(C2CCC1(C)C2(C)(C))Br",
        "C1=CC(C2C1C3(C(=C(C2(C3(Cl)Cl)Cl)Cl)Cl)Cl)Cl",
        "C=1(C(=C(C(C=1Cl)(Cl)Cl)Cl)Cl)Cl",
        "O=C1NC(=O)C(N1)(C)C",
        "C2=CC3C1C=CC(C1)C3(C2)",
        "OC(C)(CC)CC",
        "C#CC(O)(C)CC",
        "C#CC1(O)(CCCCC1)",
        "O1CCOC(C1SP(OCC)(OCC)=S)SP(OCC)(OCC)=S",
        "O=P(OCCOCCCC)(OCCOCCCC)OCCOCCCC",
        "O=C1C=C(C)CC(C)(C)C1",
        "C=CC(=C)C",
        "OCC(C)C",
        "CC(CCl)Cl",
        "NCC(N)C",
        "O=C(C)CC",
        "N#CC(O)C",
        "C(C(Cl)Cl)Cl",
        "C(=C(Cl)Cl)Cl",
        "O=C(N)C=C",
        "NNC(N)=S",
        "O=C(OC)C",
        "C(C(Cl)Cl)(Cl)Cl",
        "O=C(C=CC1=C(C)CCCC1(C)C)C",
        "Oc1c(cc(cc1Br)C(c2cc(c(O)c(c2)Br)Br)(C)C)Br",
        "Oc1c(cc(cc1Cl)C(c2cc(c(O)c(c2)Cl)Cl)(C)C)Cl",
        "Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C",
        "OC(c1ccc(cc1)Cl)(c2ccc(cc2)Cl)C",
        "Oc1ccc(cc1)C(C)(C)CC",
        "c2cc1cccc3c1c(c2)CC3",
        "c1ccc2c(c1)[nH]cc2C",
        "O=C3c5ccc1OC(C(=C)C)Cc1c5(OC4COc2cc(OC)c(OC)cc2C34)",
        "O=C(Oc1ccccc1)c3ccccc3(C(=O)Oc2ccccc2)",
        "O=C(OCC)c1ccccc1(C(=O)OCC)",
        "O=C(OCC(C)C)c1ccccc1(C(=O)OCC(C)C)",
        "O=C(OCCCC)c1ccccc1(C(=O)OCCCC)",
        "O=C(OCc1ccccc1)c2ccccc2(C(=O)OCCCC)",
        "O=[N+]([O-])c2cccc1ccccc12",
        "c1ccc3c(c1)[nH]c2ccccc23",
        "OC1OCC(O)C(O)C1(O)",
        "Oc1c(c(c(c(c1Cl)Cl)Cl)Cl)Cl",
        "O=C(OCC)C(O)C(O)C(=O)OCC",
        "Oc1c(cc(cc1Cl)Cl)Cl",
        "O=[N+]([O-])c1ccc(O)cc1C(F)(F)F",
        "O=C(N)c1ccccc1(N)",
        "O=[N+]([O-])c1ccccc1(O)",
        "O=[N+]([O-])c1cc(c(O)c(c1)C(C)CC)[N+](=O)[O-]",
        "O=[N+]([O-])c1cc(ccc1(N))C",
        "O=Cc1ccccc1(O)",
        "Oc1cccc2ccccc12",
        "O=Cc1cc(cc(c1(O))Br)Br",
        "n1cccc2ccccc12",
        "O=[N+]([O-])c1ccccc1(OC)",
        "N(CC)(CC)C1CCCCC1",
        "c1ccc(cc1)N(CC)CC",
        "OCCN(c1cccc(c1)C)CC",
        "Nc1ccc(cc1Cl)c2ccc(N)c(c2)Cl",
        "c1ccc(cc1)c2ccccc2",
        "O=C(O)C(Oc1cc(c(cc1Cl)Cl)Cl)C",
        "O=C(OCC)c1ccccc1",
        "O=C(c1ccccc1)CC(=O)C",
        "O=C(OCC)c1ccc(N)cc1",
        "ON=Cc1ccccc1(O)",
        "O=C(O)COc1ccc(cc1Cl)Cl",
        "O=C(O)CCCOc1ccc(cc1C)Cl",
        "O=C(O)CCCOc1ccc(cc1Cl)Cl",
        "O=Cc1ccc(O)cc1(O)",
        "c1ccc(c(c1)C)C",
        "Oc1ccccc1C",
        "Fc1ccccc1C",
        "Oc1ccccc1Cl",
        "c1cc(c(cc1C)C)C",
        "Oc1ccc(c(c1)C)C",
        "c1cc(c(cc1C)Cl)Cl",
        "Nc1ccc(c(c1)Cl)Cl",
        "Nc1ccc(c(N)c1)C",
        "c1c(c(cc(c1Cl)Cl)Cl)Cl",
        "Oc1cc(c(cc1Cl)Cl)Cl",
        "O=C(OCC=C)C(=C)C",
        "O1CC1c2ccccc2",
        "OCC(CBr)Br",
        "O=CC(C)CC",
        "C(C(CCl)Cl)Cl",
        "ON=C(C)CC",
        "OCCN(C(C)C)C(C)C",
        "Oc1ccc(cc1Cc2cc(ccc2(O))Cl)Cl",
        "Oc1ccc(cc1(OC))CC=C",
        "O=Cc1occc1",
        "Oc1ccc(cc1)C(C)(C)C",
        "O=C(c1ccccc1)C",
        "O=C(c1ccccc1)Cl",
        "O=[N+]([O-])c1ccccc1",
        "O=[N+]([O-])c1cccc(c1)C",
        "O=[N+]([O-])c1cc(cc(c1)[N+](=O)[O-])[N+](=O)[O-]",
        "O=[N+]([O-])c1ccc(N)c(c1)C",
        "O=Cc1cccc(c1)[N+](=O)[O-]",
        "O=[N+]([O-])c1cccc(c1)[N+](=O)[O-]",
        "c1cc(ccc1N(C)C)C",
        "O=[N+]([O-])c1ccc(cc1)C",
        "O=[N+]([O-])c1ccc(cc1)[N+](=O)[O-]",
        "c1ccc(cc1)CC",
        "C=Cc1ccccc1",
        "c1ccc(cc1)CCl",
        "c1ccc(cc1)NC",
        "ON=C1CCCCC1",
        "N#Cc1ncccc1",
        "n1ccccc1CC",
        "OCC1OC(OC1)(C)C",
        "N12(CN3CN(C1)CN(C2)C3)",
        "O=C(Nc1ccc(cc1)Cl)Nc2ccc(c(c2)Cl)Cl",
        "O(c1ccccc1)c2ccccc2",
        "c1ccc(cc1)NC(Nc2ccccc2)=S",
        "c1cc(cc(c1)C)NCC",
        "OCCN(CCO)CCO",
        "OC(C)(C)CCc1ccccc1",
        "OCCN1CCNCC1",
        "Oc1ccc(cc1)CCCCCCCCC",
        "O=C1OC(CC1)CCCCCCC",
        "OCC(CC)CCCC",
        "O=Cc1ccc(cc1)Cl",
        "n1cc(ccc1C)CC",
        "O=C(C)CCCN(CC)CC",
        "O=C(OCC)CC(=O)OCC",
        "O=C(OCCCC)C=CC(=O)OCCCC",
        "O=C(OCCCC)CCCCC(=O)OCCCC",
        "OCC=C(C)CCC=C(C)C",
        "c1cc(ccc1Br)Br",
        "Nc1ccc(cc1)Br",
        "Oc1ccc(cc1)C",
        "Oc1ccc(cc1)Cl",
        "O=C(OCC(C)C)C=C",
        "O=CC=C",
        "C=CCCl",
        "C(CCl)Cl",
        "NCCC",
        "N#CCC",
        "N#CC=C",
        "N#CCCl",
        "NCCN",
        "OCC=C",
        "C#CCO",
        "OCCO",
        "O=CC=O",
        "ON=CC",
        "OC(C)CC(O)(C)C",
        "NC(C)(C)CC(C)(C)C",
        "CC(C)(C)SC(C)(C)C",
        "O=P(OCC)(OCC)OP(=O)(OCC)OCC",
        "O=C(C)CCC",
        "O=C(OC=C)C",
        "O=C(C)CC(C)C",
        "O(C(C)C)C(C)C",
        "c1cc(cc(c1)C)C",
        "Oc1cccc(c1)C",
        "Oc1cccc(O)c1",
        "c1ccc(cc1)Br",
        "c1ccc(cc1)C",
        "n1ccc(cc1)C",
        "c1ccc(cc1)Cl",
        "O=C1CCCCC1",
        "Oc1ccccc1",
        "n1cccc(c1)C",
        "N1CCN(C)CC1",
        "n1ccccc1C",
        "N1CCNC(C)C1",
        "O=C(OCCCC)CCC",
        "O=C(O)CCCC",
        "O=C(OCCC)C",
        "C(CBr)CBr",
        "CCCCBr",
        "NCCCC",
        "N#CCC=C",
        "N#CCC#N",
        "O(C)CCN",
        "O(C)COC",
        "N(CC)CC",
        "c1cc[nH]c1",
        "O1CCCC1",
        "o1cccc1",
        "CC(C)(C)SSC(C)(C)C",
        "O=C(C)CCC(C)C",
        "O=C(OCC)CCCCCCCCC(=O)OCC",
        "O=C(C)CCCCC",
        "CCCCCC",
        "C(CCCl)CCl",
        "O=CCCCC",
        "OCCNCC",
        "C1CCCCC1",
        "O1COCOC1",
        "O=C(C)CCC=C(C)C",
        "CCCCCCBr",
        "OCCOCCO",
        "CCCSCCC",
        "OCCCCCCC",
        "CCCCCCCCBr",
        "NCCCCCCCC",
        "OCCCCCCCC",
        "O(COCCCl)CCCl",
        "O=C(C)CCCCCCCCC",
        "OCCOCCOCCO",
        "OCCCCCCCCCCC",
        "OCCCCCCCCCCCC",
        "O=C(O)CCCCCCCC=CCCCCCCCC",
        "C#CC(O)(C)C",
        "OCC(Cl)(Cl)Cl",
        "O=S1OCC2C(CO1)C3(C(=C(C2(C3(Cl)Cl)Cl)Cl)Cl)Cl",
        "OC(c1ccc(cc1)Cl)(c2ccc(cc2)Cl)C(Cl)(Cl)Cl",
        "O=S(c1ccc(OP(OCC)(OCC)=S)cc1)C",
        "O=C(ON=CC(C)(C)SC)NC",
        "O=C1c2ccccc2(C(=O)C(=C1Cl)Cl)",
        "O=C(OCC(CC)CCCC)c1ccccc1(C(=O)OCC(CC)CCCC)",
        "O=C(Oc1ccccc1)c2ccccc2(O)",
        "O=C(OCC)c1ccccc1(O)",
        "Oc1c(cc(cc1Br)Br)Br",
        "O=[N+]([O-])c1cc(N)ccc1(O)",
        "O=C(c1ccccc1)c2ccccc2",
        "OCCN(c1ccccc1)CCO",
        "O=Cc1ccc(cc1)N(CC)CC",
        "O=S(CCCCCCCC)C(C)Cc1ccc2OCOc2(c1)",
        "Oc1ccccc1(O)",
        "c1cc(c(cc1Cl)Cl)Cl",
        "O=[N+]([O-])c1ccc(c(c1)[N+](=O)[O-])C",
        "O=Cc1ccc(O)c(OCC)c1",
        "O=S(=O)(O)c1ccc(N)cc1",
        "O=[N+]([O-])c1cccc(c1)Cl",
        "O=C(OCC)CC(C(=O)OCC)SP(OC)(OC)=S",
        "O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1",
        "O=[N+]([O-])c1ccc(N)c(c1)Cl",
        "O=[N+]([O-])c1ccc(OP(OC)(OC)=S)cc1C",
        "c1ccc(cc1)Nc2ccccc2",
        "OCCOc1ccccc1",
        "O=C(OCC)CCC(=O)OCC",
        "O=C(C)CC(=O)C",
        "O=C(OCC)CCCCC",
        "O=CCCC",
        "O=C(OCCCC)C",
        "O1CCOCC1",
        "O=C(O)CCCCC(=O)O",
        "NCCCCCCCCCCCC",
        "O=[N+]([O-])C(CO)(CO)CO",
        "O=P(OCCCC)(OCCCC)OCCCC",
        "O=C1CC(=O)CC(C)(C)C1",
        "OC(C)CCl",
        "ON=C(C)C",
        "C(=C(Cl)Cl)(Cl)Cl",
        "C#CC(O)(c1ccccc1)C",
        "O=C1N(C(=O)C2CC=CCC12)SC(Cl)(Cl)Cl",
        "O=C1c2ccccc2(C(=O)N1SC(Cl)(Cl)Cl)",
        "O=C(Oc1ccccc1)c2ccc(N)cc2(O)",
        "O=C(c1cccc(c1)C)N(CC)CC",
        "Oc1ccc2ccccc2(c1)",
        "C=C(C)C1CC=C(C)CC1",
        "NCCN1CCNCC1",
        "O=C(OCC)C=C",
        "O=C(OCCCC)CCC(=O)OCCCC",
        "O=C(OCC)C=CC(=O)OCC",
        "O=C(OCC)CCCCC(=O)OCC",
        "O1C(C)CNCC1C",
        "O=C(O)CCCCC",
        "O(CCCC)CCCC",
        "O=C1C2(C3(C4(C1(C5(C2(C3(C(C45Cl)(Cl)Cl)Cl)Cl)Cl)Cl)Cl)Cl)Cl",
        "O=C(O)C2C1OC(CC1)C2(C(=O)O)",
        "c1ccc2c(c1)NC(=S)S2",
        "OCC(C)C(O)CC",
        "Oc1cccc(OC)c1",
        "Oc1ccc(OC)cc1",
        "O(c1ccc(OC)cc1)C",
        "O=P(OP(=O)(N(C)C)N(C)C)(N(C)C)N(C)C",
        "c1ccc2c(c1)c4cccc3cccc2c34",
        "n1ncc2ccccc2(c1)",
        "n1c3ccccc3(cc2ccccc12)",
        "o1ccc2ccccc12",
        "C1=CC=C2C=CC=C2(C=C1)",
        "N12(CCN(CC1)CC2)",
        "C1C2CC3CC1CC(C2)C3",
        "O=[N+]([O-])c1ccc(OP(OC)(OC)=S)cc1",
        "O(CC)P(OCC)(=S)SCSCC",
        "O(CC)P(OCC)(=S)SCCSCC",
        "O(c1cc(c(cc1Cl)Cl)Cl)P(OC)(OC)=S",
        "C1=CC2CC1C3C2C4(C(=C(C3(C4(Cl)Cl)Cl)Cl)Cl)Cl",
        "O=[N+]([O-])c1ccc(OP(=O)(OCC)OCC)cc1",
        "O=C1NC(=C(C(=O)N1C(C)CC)Br)C",
        "O=C(Oc1cc(c(c(c1)C)N(C)C)C)NC",
        "O(c1cc(c(cc1Cl)Cl)Cl)P(OCC)(CC)=S",
        "O=[N+]([O-])c1ccc(c(O)c1)[N+](=O)[O-]",
        "O=C(Nc1ccc(c(c1)Cl)Cl)N(C)C",
        "Fc2ccc(Oc1ccc(F)cc1)cc2",
        "n1c(OP(OCC)(OCC)=S)cc(nc1C(C)C)C",
        "O=[N+]([O-])c1ccc(F)cc1",
        "N#Cc1cccc(c1)C(F)(F)F",
        "Fc1ccc(N)cc1",
        "O=C(OCC)C(F)(F)F",
        "O=C(O)c1c(F)cccc1(F)",
        "O=Cc1c(F)cccc1Cl",
        "Fc1ccc(N)c(c1)C(F)(F)F",
        "N#Cc1ccccc1C(F)(F)F",
        "O=Cc1cccc(c1)C(F)(F)F",
        "OC1CC2CCC1(C)C2(C)(C)",
        "O=C1CC2CCC1(C)C2(C)(C)",
        "O1C(C)(C)C2CCC1(C)CC2",
        "O=C(O)C2(C)(CCCC3(C)(C1C(=CC(=C(C)C)CC1)CCC23))",
        "C=C2C1CCC3C1C(C)(C)CCCC23(C)",
        "c1ccc2c(c1)CCC2",
        "O2c1ccccc1CC2",
        "OC2CC1CCC2(C1)",
        "C1=CC2CCC1C2",
        "O=Cc1cnccc1",
        "O=C(CCCC)CCCC",
        "O=C(C)Br",
        "C=C(C(=C)C)C",
        "O=C1C=C(Oc2ccccc12)c3ccccc3",
        "O=[N+]([O-])c1cc(c(O)c(c1)C)[N+](=O)[O-]",
        "c1ccc(cc1)CCCCC",
        "O=C(OC(C)(C)C)C",
        "c1cc(cc(c1)Cl)Cl",
        "O=C(c1ccc(OC)cc1(O))C",
        "O=Cc1ccccc1[N+](=O)[O-]",
        "O=Cc1ccc(cc1)[N+](=O)[O-]",
        "O(CC)P(OCC)(=S)SCSP(OCC)(OCC)=S",
        "O=C(C)C(C)C",
        "O=[N+]([O-])c1cccc(c1(N))C",
        "O=[N+]([O-])c1cccc(c1(O))[N+](=O)[O-]",
        "O=[N+]([O-])c1ccc(cc1(N))C",
        "c1ccc(c(c1)Br)Br",
        "O=C=Nc1ccc(c(N=C=O)c1)C",
        "C1=C(C)CCC(=C(C)C)C1",
        "C=CCNc1ccccc1",
        "O=CCC(C)C",
        "O=C(C)CCCC",
        "C(C=CC)=CC",
        "O=C(C)CCCCCCCCCCC",
        "OC(C(C)C)C(C)C",
        "O=[N+]([O-])c1cccc(N)c1C",
        "O=[N+]([O-])c1cccc(c1C)[N+](=O)[O-]",
        "O=CN(c1ccccc1)c2ccccc2",
        "O=C(OCC)C(C(=O)OCC)Cc1ccccc1",
        "Oc1c(c(c(c(c1Br)Br)Br)Br)Br",
        "c1c(c(c(c(c1Cl)Cl)Cl)Cl)Cl",
        "Oc1c(cc(cc1I)I)I",
        "O=[N+]([O-])c1ccc(cc1[N+](=O)[O-])C",
        "O=Cc1ccc(OC)cc1(OC)",
        "O=C(Nc1ccccc1(O))C",
        "O=[N+]([O-])c1cc(c(cc1C)C)[N+](=O)[O-]",
        "O=[N+]([O-])c1cc(c(cc1(O))C)[N+](=O)[O-]",
        "O=[N+]([O-])c1cc(OCC)ccc1(N)",
        "O=[N+]([O-])c1cc(N)cc(c1)[N+](=O)[O-]",
        "O=[N+]([O-])c1ccc(c(c1)C)[N+](=O)[O-]",
        "N#Cc1cccc(c1)[N+](=O)[O-]",
        "O=C(OC)c1ccc(cc1)[N+](=O)[O-]",
        "N#Cc1ccc(cc1)[N+](=O)[O-]",
        "O=[N+]([O-])c2ccc(Oc1ccccc1)cc2",
        "O=S(Cc1ccccc1)Cc2ccccc2",
        "O=C(Nc1cccc(O)c1)C",
        "c1cc(ccc1CCl)CCl",
        "O=C(OCC)C=CC(=O)OCC",
        "o1c(ccc1C)C",
        "OCCCCl",
        "C(CCCl)CCCl",
        "CCCCCCCBr",
        "CCCSSCCC",
        "N#CCCCCCCC#N",
        "O=C(O)C(Br)Br",
        "c1cc(c(c(c1Cl)Cl)Cl)Cl",
        "Nc1ccc(c(c1Cl)Cl)Cl",
        "O=Cc1cc(ccc1(O))Cl",
        "O=Cc1c(F)c(F)c(F)c(F)c1(F)",
        "O=C(N)C(Cl)Cl",
        "NC(C)CCCCCC",
        "O=C(C)CCCCCCCC",
        "O(CCCCC)CCCCC",
        "n1cocc1C",
        "n1cc[nH]c1C",
        "O=[N+]([O-])c1ccc(cc1(O))C",
        "O=C1C2CC3CC1CC(C2)C3",
        "O=C1OC(CC1)CCCCCC",
        "O=C(Nc1ccc(c(c1)Cl)Cl)CC",
        "O=C1c2ccccc2(C(=O)N1CSP(OC)(OC)=S)",
        "Oc1c(cc(cc1C(C)(C)C)C(C)(C)C)C(C)(C)C",
        "O=CN(CCCC)CCCC",
        "NC23(CC1CC(CC(C1)C2)C3)",
        "N#CC=1C(=O)NC(=CC=1C)C",
        "O(CC)P(OCC)(=S)SCSc1ccc(cc1)Cl",
        "O=P(c1ccccc1)(c2ccccc2)c3ccccc3",
        "O=P(O)(OC)OC",
        "O=C(OCCO)C=C",
        "C#CC(O)CCCCC",
        "O=C(C)CCCCCCC",
        "C1CCC(C(C1)Cl)Cl",
        "Oc2ccc(Oc1ccccc1)cc2",
        "n1c(nc(nc1NC(C)C)SC)NCC",
        "O=C(OCCO)C(=C)C",
        "O=P(OC)OC",
        "O=Cc1ccc(cc1Cl)Cl",
        "c1ccc(cc1)SSc2ccccc2",
        "FC(F)(F)C(O)C(F)(F)F",
        "OC(C=C)CC=C",
        "C#CCCO",
        "O=C(c1cccn1C)C",
        "O=C(c1ccc(cc1)Cl)CCl",
        "n1ccc(cc1)c2ccccc2",
        "O(CC)P(CC)(=S)Sc1ccccc1",
        "O=S(c1ccccc1)c2ccccc2",
        "o1c(ccc1c2ccccc2)c3ccccc3",
        "O=C(N(C)C)C(c1ccccc1)c2ccccc2",
        "O=C(OCC(O)C)C=C",
        "O=C(O)CNCP(=O)(O)O",
        "n1cc(ccc1N)Br",
        "O=P(OCC)(OCC)Cc1ccccc1",
        "OP(OC)(OC)=S",
        "O=C(c1ccncc1)C",
        "O=C(OC)c1ccc(cc1)Cl",
        "O(c1ccccc1)CCCC",
        "N#Cc1c(cccc1Cl)Cl",
        "O=C(c1cc(OC)ccc1(OC))CBr",
        "O=P(OC(C)C)(OC(C)C)C",
        "C=Cn2c3ccccc3(c1ccccc12)",
        "O(c1cccc(N)c1)Cc2ccccc2",
        "O=C(Oc2cccc1c2(OC(C)(C)C1))NC",
        "O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])C(F)(F)F",
        "C=CCCCCCCC=C",
        "N#CCCOCCC#N",
        "O=C(O)C2(C)(CCCC3(c1ccc(cc1CCC23)C(C)C)(C))",
        "Oc1ccccc1CC=C",
        "O(C(=CCl)c1cc(ccc1Cl)Cl)P(OCC)(OCC)=S",
        "O=Cc1cc(ccc1(O))Br",
        "O=C(OC)C1C(C(=O)OC)C2(C(=C(C1(C2(Cl)Cl)Cl)Cl)Cl)Cl",
        "O=[N+]([O-])c1cc(cc(c1N(CC)CCCC)[N+](=O)[O-])C(F)(F)F",
        "C=C(CCl)CCl",
        "n1c(nc(nc1NC(C)C)Cl)NCC",
        "O=C(O)c1nc(c(c(N)c1Cl)Cl)Cl",
        "n1c(cccc1Cl)C(Cl)(Cl)Cl",
        "O=C(OCCCC)c1ccc(cc1)C(=O)OCCCC",
        "Oc2ccc(Oc1ccc(O)cc1)cc2",
        "O=C(N)c1c(cccc1Cl)Cl",
        "C#CC(O)C",
        "O=C(Oc1ccc(c(c1)C)N(C)C)NC",
        "n1c([nH]c(c1Br)Br)Br",
        "c1ccc(cc1)c2ccccc2Cl",
        "c1ccc(cc1)c2cccc(c2)Cl",
        "O=[N+]([O-])c2ccc(OP(OCC)(c1ccccc1)=S)cc2",
        "C#CCC(O)C",
        "Oc1ccc(cc1(O))Cl",
        "O=C(OC)c1ccc(O)cc1(O)",
        "O=C(NC1CC2CC1C3CCCC23)N(C)C",
        "n1c(c(c(c(c1Cl)Cl)Cl)Cl)Cl",
        "OC1CC(C)CCC1(C(C)C)",
        "O=S(=O)(c1ccc(cc1)C)n2cncc2",
        "O=C(c1ccc(cc1Cl)Cl)C",
        "N#CCCCCCCCC",
        "c2cc1cc(c(cc1c(c2)C)C)C",
        "Fc1ccc(N)cc1C(F)(F)F",
        "OC2CCCCC2(c1ccccc1)",
        "O=C(OCCOCC)C(=C)C",
        "C13(C4(C2(C5(C(C1(C2(Cl)Cl)Cl)(C3(C(C45Cl)(Cl)Cl)Cl)Cl)Cl)Cl)Cl)Cl",
        "Oc1c(ccc(c1C)C)C",
        "O=C(N)c1ccccc1(OC)",
        "O=C(N)c1ccc(cc1Cl)Cl",
        "O=C(OCC1OCCC1)C(=C)C",
        "Oc1cc(c(cc1(OC))Cl)Cl",
        "O=C(OCc1ccccc1)C(=C)C",
        "O=C(OCCCCCC)C=C",
        "O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C1",
        "c1cc[n+]3c(c1)c2cccc[n+]2CC3",
        "OCCCc1cnccc1",
        "NCCCCCCCCCCCCC",
        "O=C(c1ccc(cc1)Cl)c2ccccc2(N)",
        "O=C(OC)c1cc(ccc1Cl)Cl",
        "n1c(OP(OCC)(OCC)=S)c(cc(c1Cl)Cl)Cl",
        "O=Cc1cc(OC)c(O)c(c1)Br",
        "O=C(OC1CCCCC1)C=C",
        "N#CC(Br)Br",
        "O(c1ccc(cc1)Sc2ccc(OP(OC)(OC)=S)cc2)P(OC)(OC)=S",
        "Oc1cc(c(cc1(O))Cl)Cl",
        "n1c(cccc1c2ccccc2)c3ccccc3",
        "O=S(=O)(Oc1ccc(cc1Cl)Cl)C",
        "O=[N+]([O-])c1cc(c(cc1Cl)Cl)[N+](=O)[O-]",
        "C#CC(O)(c1ccccc1)c2ccccc2",
        "O=CC(C)C(C)CC",
        "OCC#CCCCCCCC",
        "O=S(=O)(NCC)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
        "O=CC=CC",
        "O(c1ccc(cc1)C=CC)C",
        "CC(C)SSC(C)C",
        "OCc1cocc1",
        "O=Cc1cc(OC)c(OC)cc1(OC)",
        "O=C(O)CC(C)(C)CC(=O)O",
        "Oc1cc(c(c(c1Cl)Cl)Cl)Cl",
        "n1ccc(cc1)CCc2ccncc2",
        "O=[N+]([O-])c1c(O)cccc1C",
        "O=C1N(C(N(C(=O)C1)CC)=S)CC",
        "c1cc2c(cc1Cl)NC(=S)S2",
        "O=C(OC)c1ccc(C(=O)OC)c(N)c1",
        "CCSCCSCC",
        "O=C(OCCCCCCCC)C(O)C",
        "O=C(c1ccc(c(c1)[N+](=O)[O-])Cl)C",
        "O=P(Oc1nc(c(cc1Cl)Cl)Cl)(OC)OC",
        "n1c(nc(cc1C)Cl)N",
        "O(c1cccc(OC)c1C)C",
        "n1ccccc1N(C)C",
        "NCC(C)(C)C",
        "O=C(O)C1(C)(CCCC2(C)(C3C(=CCC12)CC(C=C)(C)CC3))",
        "N#Cc1cc(ccc1(N))Cl",
        "C=C(C)C1CC=C(C)CC1",
        "O=C(C)CCCCCCCCCC",
        "O=CC=Cc1ccc(cc1)N(C)C",
        "O=[N+]([O-])c1c(cc(c(c1Cl)[N+](=O)[O-])Cl)Cl",
        "O=Cc1cc(ccc1Cl)[N+](=O)[O-]",
        "O=C(Nc1cccc(c1)Cl)c2cccc(c2(O))[N+](=O)[O-]",
        "N#Cc1c(cccc1Cl)C",
        "Oc1cccnc1Br",
        "Oc1cccnc1Cl",
        "C#CCN(CC#C)CC#C",
        "O(CC)C(OCC)CN(C)CC(OCC)OCC",
        "NCCCN1CCN(CC1)CCCN",
        "OC(C=C)(C)CCC=C(C)CCC=C(C)C",
        "C1=C(C)C2CC(C1)C2(C)(C)",
        "O=Cc1ccc(OCC)cc1",
        "O=[N+]([O-])c1ccc(N)c(c1C)[N+](=O)[O-]",
        "N#CC(C(=O)N)(Br)Br",
        "O(CC)P(OCC)(=S)SCSC(C)(C)C",
        "O=C(C(=C(OP(=O)(OC)OC)C)Cl)N(CC)CC",
        "c1ccc(c(c1)C(Br)Br)C(Br)Br",
        "O=[N+]([O-])c2oc(C=Cc1nc(ccc1)CO)cc2",
        "O=C(c1ccc(c(c1Cl)Cl)Cl)C",
        "O=C(c1ccc(OC)c(OC)c1(OC))C",
        "NC(C)CC",
        "c1ccc(cc1)CNCC",
        "O=C(c1cnccc1)c2ccc(cc2)Br",
        "O=C(c1ccncc1)c2ccccc2",
        "O1C(C)(C)CCC1(C)C",
        "O=C(OCC(OCCCC)C)C(Oc1cc(c(cc1Cl)Cl)Cl)C",
        "O=[N+]([O-])c1ncccc1(O)",
        "NCCC(C)(C)C",
        "O=C1N(C(=O)C(N1Br)(C)C)Cl",
        "Nc1ccc(cc1)CCCCCCCC",
        "O=P(O)(O)CCCl",
        "O=C(ON=C(C)SC)NC",
        "Oc1nc(ccc1)Cl",
        "O=C(Nc1ccc(cc1Cl)[N+](=O)[O-])c2c(O)c(cc(c2C)Cl)C(C)(C)C",
        "n1nc(c(nc1N)C)C",
        "O=Cc1ccc(cc1(O))N(CC)CC",
        "C=C1CCC2CC1C2(C)(C)",
        "O=[N+]([O-])c1ccc(c(c1[N+](=O)[O-])C)[N+](=O)[O-]",
        "n1c(cccc1Cl)C",
        "O=C(OCCOCCCC)C(Oc1cc(c(cc1Cl)Cl)Cl)C",
        "O=[N+]([O-])c1cc(N)cc(c1C)[N+](=O)[O-]",
        "C#CC(O)(C)CCC(C)C",
        "OCCN(CC=C(C)CCC=C(C)C)CC=C(C)CCC=C(C)C",
        "n1c(oc(c1C)C)C",
        "O=C(C)C(C)CN(C)C",
        "O=C(OCC)C(N(C(=O)c1ccccc1)c2ccc(c(c2)Cl)Cl)C",
        "O=C(OCc1coc(c1)Cc2ccccc2)C4C(C=C3CCCC3)C4(C)(C)",
        "O=C(N)c1cccc(c1)Br",
        "n1c(nc(nc1NC2CC2)Cl)NC(C)C",
        "O=C(OCC)CC(C(=O)OCC)S",
        "O=C(ON=C(C(=O)N(C)C)SC)NC",
        "Nc1c(cccc1C(C)C)C(C)C",
        "O=C(Oc1cccc(c1)C(C)CC)N(C)Sc2ccccc2",
        "O=C1C=CSN1CCCCCCCC",
        "O1C2C3C(C(C12Cl)Cl)C4(C(=C(C3(C4(Cl)Cl)Cl)Cl)Cl)Cl",
        "O=C(OC1C(=C(C(=O)C1)CC=C)C)C2C(C=C(C)C)C2(C)(C)",
        "O=[N+]([O-])c1cc(c(N)c(c1N(CC)CC)[N+](=O)[O-])C(F)(F)F",
        "n1c(OP(OC)(OC)=S)cc(nc1N(CC)CC)C",
        "FC(F)C(F)(F)C(O)(C)C",
        "O=C(OCCCCCC(C)C)C(Oc1cc(c(cc1Cl)Cl)Cl)C",
        "n1c(nc(nc1NCC)NC(C)(C)C)OC",
        "O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])C(C)C",
        "O=C(OCC)C(N(C(=O)c1ccccc1)c2ccc(c(c2)Cl)Cl)C",
        "O(c1cc(cc(OC)c1(OC))CNCCCOC)C",
        "O1CCCCC1CBr",
        "O=[N+]([O-])c1cc(N)c(c(c1)[N+](=O)[O-])C",
        "Nc1ccc(cc1)CCCCCCCCCC",
        "c1ccc(c(c1)c2cc(ccc2Cl)Cl)Cl",
        "O=P(Oc1ccc(cc1Cl)Br)(OCC)SCCC",
        "O=C(OC)c1ccc(cc1[N+](=O)[O-])Cl",
        "O=Cc1cc(O)ccc1[N+](=O)[O-]",
        "O=C(N(c1c(cccc1CC)C)C(C)COC)CCl",
        "N#CC(OC(=O)C(c1ccc(cc1)Cl)C(C)C)c3cccc(Oc2ccccc2)c3",
        "CCSCCCCSCC",
        "O=C(N)c1ccc(cc1)C(C)(C)C",
        "O=[N+]([O-])c1ccc(c(c1(N))C)[N+](=O)[O-]",
        "O1CC1(c2cc(cc(c2)Cl)Cl)CC(Cl)(Cl)Cl",
        "C#CC(O)C",
        "N#CC(OC(=O)C(c1ccc(cc1)Cl)C(C)C)c3cccc(Oc2ccccc2)c3",
        "O=C(Nc1cccc(OC(C)C)c1)c2ccccc2C(F)(F)F",
        "O=C(OCC)C(Oc3ccc(Oc1nc2ccc(cc2(o1))Cl)cc3)C",
        "N#CC(OC(=O)C1C(C=C(Cl)Cl)C1(C)(C))c3ccc(F)c(Oc2ccccc2)c3",
        "O=S(=O)([O-])c1ccc[n+](c1)Cc2ccccc2",
        "O=Cc2cccc(Oc1ccc(cc1)C(C)(C)C)c2",
        "N#CC(OC(=O)C(c1ccc(OC(F)F)cc1)C(C)C)c3cccc(Oc2ccccc2)c3",
        "O=[N+]([O-])c1ccc(c(c1(N))[N+](=O)[O-])C",
        "O=CNc1cccc(c1C)Cl",
        "O=Cc2cccc(Oc1ccc(c(c1)Cl)Cl)c2",
        "O=C1C(=C(C=NN1C(C)(C)C)SCc2ccc(cc2)C(C)(C)C)Cl"
    };
}
