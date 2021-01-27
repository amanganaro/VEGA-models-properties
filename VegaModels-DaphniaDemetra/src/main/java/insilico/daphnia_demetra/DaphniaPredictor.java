package insilico.daphnia_demetra;

public class DaphniaPredictor extends Predictor{


    double minTrainingSetNEGRI = -0.983789,
            maxTrainingSetNEGRI = 8.88255,
            minTrainingSetCSL_03 = -0.248085484,
            maxTrainingSetCSL_03 = 7.253535748,
            minTrainingSetCSL_04 = -0.29234001,
            maxTrainingSetCSL_04 = 7.203043938;

    /** Creates a new instance of DaphniaPredictor */
    public DaphniaPredictor(){
        min = -0.95;
        max = 7.63;
        measureUnit = "[-LOG(mmol/l)]";
        secondMeasureUnit = "[mg/l]";
/*      String[] ajut1 = {
                "MW (Dragon)",
                "Mp (Dragon)",
                "nP (Dragon)",
                "nNR2Ph (Dragon)",
                "nHAcc (Dragon)",
                "O-057 (Dragon)",
                "O-060 (Dragon)",
                "S-107 (Dragon)",
                "WA (Dragon)",
                "IDE (Dragon)",
                "IC2 (Dragon)",
                "Eig1p (Dragon)",
                "T(F..Cl) (Dragon)",
                "SRW05 (Dragon)",
                "BEHm1 (Dragon)",
                "MLOGP (Dragon)"
                         };
 */
        String[] ajut1 = {
                "BEHm1 (Dragon)",
                "Eig1p (Dragon)",
                "IC2 (Dragon)",
                "IDE (Dragon)",
                "MLOGP (Dragon)",
                "Mp (Dragon)",
                "MW (Dragon)",
                "nHAcc (Dragon)",
                "nNR2Ph (Dragon)",
                "nP (Dragon)",
                "O-057 (Dragon)",
                "O-060 (Dragon)",
                "S-107 (Dragon)",
                "SRW05 (Dragon)",
                "T(F..Cl) (Dragon)",
                "WA (Dragon)"
        };
        descrName = ajut1;
/*            String[] ajut2 = {
                "Symbol: MW; Definition: molecular weight; Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: Mp; Definition: mean atomic polarizability (scaled on Carbon atom); Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: nP; Definition: number of Phosphorous atoms; Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: nNR2Ph; Definition: number of tertiary amines (aromatic); Class: functional groups; Dimensionality: 2D",
                "Symbol: nHAcc; Definition: number of acceptor atoms for H-bonds (N O F); Class: functional groups; Dimensionality: 2D",
                "Symbol: O-057; Definition: phenol / enol / carboxyl OH; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: O-060; Definition: Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: S-107; Definition: R2S / RS-SR; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: WA; Description: mean Wiener index; Clsss: topological descriptors; Dimensionality: 2D",
                "Symbol: IDE; Definition: mean information content on the distance equality; Class: topological descriptors; Dimensionality: ",
                "Symbol: IC2; Definition: information content index (neighborhood symmetry of 2-order); Class: topological descriptors; Dimensionality: 2D",
                "Symbol: Eig1p; Definition: Leading eigenvalue from polarizability weighted distance matrix; Class: topological descriptors; Dimensionality: 2D",
                "Symbol: T(F..Cl); Definition: sum of topological distances between F..Cl; Class: topological descriptors; Dimensionality: 2D",
                "Symbol: SRW05; Definition: self-returning walk count of order 05; Class: molecular walk counts; Dimensionality: 2D",
                "Symbol: BEHm1; Definition: highest eigenvalue n. 1 of Burden matrix / weighted by atomic masses; Class: BCUT descriptors; Dimensionality: 2D",
                "Symbol: MLOGP; Definition: Moriguchi octanol-water partition coeff. (logP); Class: properties; Dimensionality: 2D"
            };
 */
        String[] ajut2 = {
                "Symbol: BEHm1; Definition: highest eigenvalue n. 1 of Burden matrix / weighted by atomic masses; Class: BCUT descriptors; Dimensionality: 2D",
                "Symbol: Eig1p; Definition: Leading eigenvalue from polarizability weighted distance matrix; Class: topological descriptors; Dimensionality: 2D",
                "Symbol: IC2; Definition: information content index (neighborhood symmetry of 2-order); Class: topological descriptors; Dimensionality: 2D",
                "Symbol: IDE; Definition: mean information content on the distance equality; Class: topological descriptors; Dimensionality: ",
                "Symbol: MLOGP; Definition: Moriguchi octanol-water partition coeff. (logP); Class: properties; Dimensionality: 2D",
                "Symbol: Mp; Definition: mean atomic polarizability (scaled on Carbon atom); Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: MW; Definition: molecular weight; Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: nHAcc; Definition: number of acceptor atoms for H-bonds (N O F); Class: functional groups; Dimensionality: 2D",
                "Symbol: nNR2Ph; Definition: number of tertiary amines (aromatic); Class: functional groups; Dimensionality: 2D",
                "Symbol: nP; Definition: number of Phosphorous atoms; Class: constitutional descriptors; Dimensionality: 2D",
                "Symbol: O-057; Definition: phenol / enol / carboxyl OH; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: O-060; Definition: Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: S-107; Definition: R2S / RS-SR; Class: atom-centred fragments; Dimensionality: 2D",
                "Symbol: SRW05; Definition: self-returning walk count of order 05; Class: molecular walk counts; Dimensionality: 2D",
                "Symbol: T(F..Cl); Definition: sum of topological distances between F..Cl; Class: topological descriptors; Dimensionality: 2D",
                "Symbol: WA; Description: mean Wiener index; Clsss: topological descriptors; Dimensionality: 2D"
        };
        descrDetails = ajut2;
    }

    public double maximumSubmodelsPrediction(double[] input) {
        //int[] negriDescriptors = {1,2,4,5,6,9,10,12,13,14,15}; //old index
        int[] negriDescriptors = {5,9,7,10,11,3,2,14,13,0,4}; //new sorted index
        double[] inputNegri = new double[negriDescriptors.length];
        for(int i=0;i<negriDescriptors.length;i++)
            inputNegri[i] = input[negriDescriptors[i]];
        double negri = negri_d05(inputNegri);

        //int[] csl_d03Descriptors = {0,1,2,3,5,6,7,11,13,14,15}; // old index
        int[] csl_d03Descriptors = {6,5,9,8,10,11,12,1,13,0,4}; // new sorted index
        double[] inputCsl_d03 = new double[csl_d03Descriptors.length];
        for(int i=0;i<csl_d03Descriptors.length;i++)
            inputCsl_d03[i] = input[csl_d03Descriptors[i]];
        double csl_d03 = csl_d03(inputCsl_d03);

        //int[] csl_d04Descriptors = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14}; // old index
        int[] csl_d04Descriptors = {6,5,9,8,7,10,11,12,15,3,2,1,14,13,0};// new sorted index
        double[] inputCsl_d04 = new double[csl_d04Descriptors.length];
        for(int i=0;i<csl_d04Descriptors.length;i++)
            inputCsl_d04[i] = input[csl_d04Descriptors[i]];
        double csl_d04 = csl_d04(inputCsl_d04);

        double valueNEGRI = 0, valueCSL_03 = 0, valueCSL_04 = 0;

        // NEGRI
        if(negri < minTrainingSetNEGRI)
            valueNEGRI = minTrainingSetNEGRI;
        else if(negri > maxTrainingSetNEGRI)
            valueNEGRI = maxTrainingSetNEGRI;
        else
            valueNEGRI = negri;

        // CSL_03
        if(csl_d03 < minTrainingSetCSL_03)
            valueCSL_03 = minTrainingSetCSL_03;
        else if(csl_d03 > maxTrainingSetCSL_03)
            valueCSL_03 = maxTrainingSetCSL_03;
        else
            valueCSL_03 = csl_d03;

        // CSL_04
        if(csl_d04 < minTrainingSetCSL_04)
            valueCSL_04 = minTrainingSetCSL_04;
        else if(csl_d04 > maxTrainingSetCSL_04)
            valueCSL_04 = maxTrainingSetCSL_04;
        else
            valueCSL_04 = csl_d04;

        return Math.max(valueNEGRI,Math.max(valueCSL_03, valueCSL_04));
    }

    public double minimumSubmodelsPrediction(double[] input) {
        //int[] negriDescriptors = {1,2,4,5,6,9,10,12,13,14,15}; //old index
        int[] negriDescriptors = {5,9,7,10,11,3,2,14,13,0,4}; //new sorted index
        double[] inputNegri = new double[negriDescriptors.length];
        for(int i=0;i<negriDescriptors.length;i++)
            inputNegri[i] = input[negriDescriptors[i]];
        double negri = negri_d05(inputNegri);

        //int[] csl_d03Descriptors = {0,1,2,3,5,6,7,11,13,14,15}; // old index
        int[] csl_d03Descriptors = {6,5,9,8,10,11,12,1,13,0,4}; // new sorted index
        double[] inputCsl_d03 = new double[csl_d03Descriptors.length];
        for(int i=0;i<csl_d03Descriptors.length;i++)
            inputCsl_d03[i] = input[csl_d03Descriptors[i]];
        double csl_d03 = csl_d03(inputCsl_d03);

        //int[] csl_d04Descriptors = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14}; // old index
        int[] csl_d04Descriptors = {6,5,9,8,7,10,11,12,15,3,2,1,14,13,0};// new sorted index
        double[] inputCsl_d04 = new double[csl_d04Descriptors.length];
        for(int i=0;i<csl_d04Descriptors.length;i++)
            inputCsl_d04[i] = input[csl_d04Descriptors[i]];
        double csl_d04 = csl_d04(inputCsl_d04);

        double valueNEGRI = 0, valueCSL_03 = 0, valueCSL_04 = 0;

        // NEGRI
        if(negri < minTrainingSetNEGRI)
            valueNEGRI = minTrainingSetNEGRI;
        else if(negri > maxTrainingSetNEGRI)
            valueNEGRI = maxTrainingSetNEGRI;
        else
            valueNEGRI = negri;

        // CSL_03
        if(csl_d03 < minTrainingSetCSL_03)
            valueCSL_03 = minTrainingSetCSL_03;
        else if(csl_d03 > maxTrainingSetCSL_03)
            valueCSL_03 = maxTrainingSetCSL_03;
        else
            valueCSL_03 = csl_d03;

        // CSL_04
        if(csl_d04 < minTrainingSetCSL_04)
            valueCSL_04 = minTrainingSetCSL_04;
        else if(csl_d04 > maxTrainingSetCSL_04)
            valueCSL_04 = maxTrainingSetCSL_04;
        else
            valueCSL_04 = csl_d04;

        return Math.min(valueNEGRI,Math.min(valueCSL_03, valueCSL_04 ));
    }

    public double predict(double[] input) {
        double result = 0;
        //int[] negriDescriptors = {1,2,4,5,6,9,10,12,13,14,15}; //old index
        int[] negriDescriptors = {5,9,7,10,11,3,2,14,13,0,4}; //new sorted index
        double[] inputNegri = new double[negriDescriptors.length];
        for(int i=0;i<negriDescriptors.length;i++)
            inputNegri[i] = input[negriDescriptors[i]];
        double negri = negri_d05(inputNegri);

        //int[] csl_d03Descriptors = {0,1,2,3,5,6,7,11,13,14,15}; // old index
        int[] csl_d03Descriptors = {6,5,9,8,10,11,12,1,13,0,4}; // new sorted index
        double[] inputCsl_d03 = new double[csl_d03Descriptors.length];
        for(int i=0;i<csl_d03Descriptors.length;i++)
            inputCsl_d03[i] = input[csl_d03Descriptors[i]];
        double csl_d03 = csl_d03(inputCsl_d03);

        //int[] csl_d04Descriptors = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14}; // old index
        int[] csl_d04Descriptors = {6,5,9,8,7,10,11,12,15,3,2,1,14,13,0};// new sorted index
        double[] inputCsl_d04 = new double[csl_d04Descriptors.length];
        for(int i=0;i<csl_d04Descriptors.length;i++)
            inputCsl_d04[i] = input[csl_d04Descriptors[i]];
        double csl_d04 = csl_d04(inputCsl_d04);

        double average = ((negri+csl_d03+csl_d04)/3),
                minimum = Math.min(negri,Math.min(csl_d03, csl_d04)),
                maximum= Math.max(negri,Math.max(csl_d03, csl_d04));
//System.out.println(average+"  "+"  "+minimum+"  "+maximum);
        if (average>6.20698945995121)
            result =1.05332910664321*minimum+0.123719191493296;
        else if (average>1.85482544979933)
            result=1.0642232289387*average+0.0522008872240036;
        else
            result = 0.968589465067067*average-0.0351616731885144;
        return result;
    }

    private double negri_d05(double[] input){
        double[] negriCoeff = {2.35041,1.19186,-0.0353175,-0.929535,0.599448,0.53352,-0.470345,0.0239859,0.0241303,0.264749,0.20342};
        double result = -1.04056;
        for(int i=0;i<negriCoeff.length;i++)
            result+=negriCoeff[i]*input[i];
        return result;
    }

    private double csl_d03(double[] input){
        double[][] cslNormCoeff={
                {0.001837353222072,	-1.034365653991690},
                {2.950819730758660,	-2.493442773818970},
                {0.899999976158142,	-0.899999976158142},
                {1.799999952316280,	-0.899999976158142},
                {0.599999964237213,	-0.899999976158142},
                {0.299999982118606,	-0.899999976158142},
                {0.599999964237213,	-0.899999976158142},
                {0.004469673149288,	-0.903691947460174},
                {0.022499999031425,	-0.899999976158142},
                {0.256629586219787,	-1.777929782867430},
                {0.137551575899124,	-0.611554324626922}
        };
        double[][] cslhiddLayerWeigths={
                {0.424578636884689,	-0.259880155324935,	1.114727258682250,	1.258683681488030,	0.361701726913452},
                {-0.570167005062103,	-0.044166374951601,	-0.312235087156295,	-0.777282953262329,	0.112529709935188},
                {-0.290090203285217,	-0.295901954174041,	-0.421184360980987,	-0.028796726837754,	-0.955345153808593},
                {-0.066141337156296,	-0.344934463500976,	0.374694228172302,	-0.234412059187889,	-0.458430081605911},
                {0.395887076854705,	0.458881586790084,	-0.253779441118240,	0.716379463672637,	1.208142995834350},
                {-1.081686139106750,	0.268210858106613,	-0.195189177989959,	-0.646822750568389,	-0.559111952781677},
                {0.343890577554702,	0.133807271718978,	-0.020534679293633,	-0.361936241388320,	0.331734657287597},
                {1.242887020111080,	-0.332801312208175,	0.463768512010574,	-0.750640571117401,	-0.195066854357719},
                {0.237126320600509,	0.108352057635784,	0.389396011829376,	-0.095403492450714,	-0.625398933887481},
                {0.367807388305664,	0.404628574848175,	-0.541203558444976,	0.440541654825210,	-1.373045206069940},
                {1.327612280845640,	0.043651398271322,	0.355900347232818,	-0.780174791812896,	-0.037750490009785},
                {0.169228896498680,	-0.191277042031288,	0.293937206268310,	-0.737077713012695,	-0.267524421215057}
        };
        double[] cslOutLayerWeigths={0.392974376678466,0.069456756114960,-0.286341369152069,-0.694058954715728,-0.589089572429656,0.309618264436721};
        double[] denormWeigths={0.209796354174614,	-0.700214982032775};

        double result = 0;
        double[] normInput = new double[input.length];
        for(int i=0;i<input.length;i++)
            normInput[i]=input[i]*cslNormCoeff[i][0]+cslNormCoeff[i][1];
        double[] hiddLayerOutput = new double[cslhiddLayerWeigths[0].length];
        for(int i=0;i<cslhiddLayerWeigths[0].length;i++){
            hiddLayerOutput[i] = cslhiddLayerWeigths[cslhiddLayerWeigths.length-1][i];
            for(int j=0;j<normInput.length;j++)
                hiddLayerOutput[i]+=normInput[j]*cslhiddLayerWeigths[j][i];
        }
        for(int i =0;i<hiddLayerOutput.length;i++)
            hiddLayerOutput[i]=(Math.exp(hiddLayerOutput[i])-Math.exp(-hiddLayerOutput[i]))/(Math.exp(hiddLayerOutput[i])+Math.exp(-hiddLayerOutput[i]));
        double outLayerOutput=cslOutLayerWeigths[cslOutLayerWeigths.length-1];
        for(int i=0;i<hiddLayerOutput.length;i++)
            outLayerOutput+=hiddLayerOutput[i]*cslOutLayerWeigths[i];
        outLayerOutput=(Math.exp(outLayerOutput)-Math.exp(-outLayerOutput))/(Math.exp(outLayerOutput)+Math.exp(-outLayerOutput));
        result = 1/denormWeigths[0]*outLayerOutput-(denormWeigths[1]/denormWeigths[0]);
        return result;
    }
    private double csl_d04(double[] input){
        double[][] cslNormCoeff={
                {0.001837353222072,	-1.034365653991690},
                {2.950819730758660,	-2.493442773818970},
                {0.899999976158142,	-0.899999976158142},
                {1.799999952316280,	-0.899999976158142},
                {0.089999996125698,	-0.899999976158142},
                {0.599999964237213,	-0.899999976158142},
                {0.299999982118606,	-0.899999976158142},
                {0.599999964237213,	-0.899999976158142},
                {0.215853214263916,	-1.115853190422050},
                {0.407885760068893,	-0.899999976158142},
                {0.566393971443176,	-1.676526188850400},
                {0.004469673149288,	-0.903691947460174},
                {0.021428570151329,	-0.899999976158142},
                {0.022499999031425,	-0.899999976158142},
                {0.256629586219787,	-1.777929782867430}
        };
        double[][] cslhiddLayerWeigths={
                {0.537888944149017,	0.372500747442245,	0.302290856838226,	-0.413798272609710,	0.276089340448379},
                {0.388986080884933,	0.055265896022320,	-1.413969874382010,	-0.400761008262634,	-0.135487109422683},
                {-0.667276740074157,	0.016729231923819,	-0.273116886615753,	0.033737663179636,	0.055062416940928},
                {0.342379897832870,	0.245219349861145,	-0.597729027271270,	0.422525405883789,	0.047544516623020},
                {-0.008485855534673,	-0.183579742908477,	0.757963955402374,	0.125155672430992,	0.905733525753021},
                {0.415688008069992,	0.138378307223320,	0.761967599391937,	-0.651618719100952,	-0.567580342292785},
                {-0.403623729944229,	-0.057683683931828,	-0.803381919860839,	0.754199743270874,	0.629846692085266},
                {0.089607253670692,	0.340981423854827,	0.546277940273284,	0.144883304834365,	-0.766449809074401},
                {0.319028079509735,	0.550674140453338,	1.071244120597830,	0.652613639831542,	-0.266038447618484},
                {-0.621404588222503,	-0.387048453092575,	0.303762733936309,	1.154203772544860,	0.110996283590793},
                {0.484364897012710,	-1.262500166893000,	0.496947497129440,	-0.242692589759826,	-0.275901913642883},
                {-0.594185054302215,	-0.076654158532619,	0.235372602939605,	1.184027433395380,	-0.807588875293731},
                {-0.593708097934722,	0.478569209575653,	0.191449463367462,	0.377512305974960,	-0.418107718229293},
                {-0.592061161994934,	0.091559417545795,	-0.058782465755939,	-0.514451146125793,	-0.586326479911804},
                {-1.164556622505180,	-0.351711362600326,	-0.670713603496551,	-0.087889552116394,	0.949004769325256},
                {-0.150861650705337,	-0.259994089603424,	-0.156676724553108,	0.824186563491821,	-0.099451787769794}
        };
        double[] cslOutLayerWeigths={-0.500746428966522,-0.602665007114410,-0.641864597797393,0.571030974388122,-0.541331589221954,-0.030191071331501};
        double[] denormWeigths={0.209796354174614,-0.700214982032775};

        double result = 0;
        double[] normInput = new double[input.length];
        for(int i=0;i<input.length;i++)
            normInput[i]=input[i]*cslNormCoeff[i][0]+cslNormCoeff[i][1];
        double[] hiddLayerOutput = new double[cslhiddLayerWeigths[0].length];
        for(int i=0;i<cslhiddLayerWeigths[0].length;i++){
            hiddLayerOutput[i] = cslhiddLayerWeigths[cslhiddLayerWeigths.length-1][i];
            for(int j=0;j<normInput.length;j++)
                hiddLayerOutput[i]+=normInput[j]*cslhiddLayerWeigths[j][i];
        }
        for(int i =0;i<hiddLayerOutput.length;i++)
            hiddLayerOutput[i]=(Math.exp(hiddLayerOutput[i])-Math.exp(-hiddLayerOutput[i]))/(Math.exp(hiddLayerOutput[i])+Math.exp(-hiddLayerOutput[i]));
        double outLayerOutput=cslOutLayerWeigths[cslOutLayerWeigths.length-1];
        for(int i=0;i<hiddLayerOutput.length;i++)
            outLayerOutput+=hiddLayerOutput[i]*cslOutLayerWeigths[i];
        outLayerOutput=(Math.exp(outLayerOutput)-Math.exp(-outLayerOutput))/(Math.exp(outLayerOutput)+Math.exp(-outLayerOutput));
        result = 1/denormWeigths[0]*outLayerOutput-(denormWeigths[1]/denormWeigths[0]);
        return result;
    }
}
