package insilico.skin_irfmn.desciptors.weights;

import insilico.core.descriptor.Descriptor;
import insilico.core.descriptor.DescriptorBlock;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.tools.utils.MoleculeUtilities;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.ArrayList;
import java.util.List;

public class FG extends DescriptorBlock {

    private static final long serialVersionUID = 1L;
    private static final String BlockName = "Functional Groups";

//    private QueryAtomContainer[] Queries; // one for each pre-fetched SMARTS
    private Pattern[] Queries;


    // Definition of FG names, description and SMARTS

    private static String[] FGNames = {
            "nCp",
            "nCs",
            "nCt",
            "nCq",
            "nCrs",
            "nCrt",
            "nCrq",
            "nCar",
            "nCbH",
            "nCb–",
            "nCconj",
            "nR=Cp",
            "nR=Cs",
            "nR=Ct",
            "n=C=",
            "nR#CH/X",
            "nR#C–",
            "nROCN",
            "nArOCN",
            "nRNCO",
            "nArNCO",
            "nRSCN",
            "nArSCN",
            "nRNCS",
            "nArNCS",
            "nRCOOH",
            "nArCOOH",
            "nRCOOR",
            "nArCOOR",
            "nRCONH2",
            "nArCONH2",
            "nRCONHR",
            "nArCONHR",
            "nRCONR2",
            "nArCONR2",
            "nROCON",
            "nArOCON",
            "nRCOX",
            "nArCOX",
            "nRCSOH",
            "nArCSOH",
            "nRCSSH",
            "nArCSSH",
            "nRCOSR",
            "nArCOSR",
            "nRCSSR",
            "nArCSSR",
            "nRCHO",
            "nArCHO",
            "nRCO",
            "nArCO",
            "nCONN",
            "nC=O(OR)2",
            "nN=C-N<",
            "nC(=N)N2",
            "nRC=N",
            "nArC=N",
            "nRCNO",
            "nArCNO",
            "nRNH2",
            "nArNH2",
            "nRNHR",
            "nArNHR",
            "nRNR2",
            "nArNR2",
            "nN-N",
            "nN=N",
            "nRCN",
            "nArCN",
            "nN+",
            "nNq",
            "nRNHO",
            "nArNHO",
            "nRNNOx",
            "nArNNOx",
            "nRNO",
            "nArNO",
            "nRNO2",
            "nArNO2",
            "nN(CO)2",
            "nC=N-N<",
            "nROH",
            "nArOH",
            "nOHp",
            "nOHs",
            "nOHt",
            "nROR",
            "nArOR",
            "nROX",
            "nArOX",
            "nO(C=O)2",
            "nH2O",
            "nSH",
            "nC=S",
            "nRSR",
            "nRSSR",
            "nSO",
            "nS(=O)2",
            "nSOH",
            "nSOOH",
            "nSO2OH",
            "nSO3OH",
            "nSO2",
            "nSO3",
            "nSO4",
            "nSO2N",
            "nPO3",
            "nPO4",
            "nPR3",
            "nP(=O)O2R",
            "nP(=O)R3/nPR5",
            "nCH2RX",
            "nCHR2X",
            "nCR3X",
            "nR=CHX",
            "nR=CRX",
            "nR#CX",
            "nCHRX2",
            "nCR2X2",
            "nR=CX2",
            "nCRX3",
            "nArX",
            "nCXr",
            "nCXr=",
            "nCconjX",
            "nAziridines",
            "nOxiranes",
            "nThiranes",
            "nAzetidines",
            "nOxetanes",
            "nThioethanes",
            "nBeta-Lactams",
            "nPyrrolidines",
            "nOxolanes",
            "nth-Thiophenes",
            "nPyrroles",
            "nPyrazoles",
            "nImidazoles",
            "nFuranes",
            "nThiophenes",
            "nOxazoles",
            "nIsoxazoles",
            "nThiazoles",
            "nIsothiazoles",
            "nTriazoles",
            "nPyridines",
            "nPyridazines",
            "nPyrimidines",
            "nPyrazines",
            "n135-Triazines",
            "n124-Triazines",
            "nHDon",
            "nHAcc",
            "nHBonds"
    };

    private static String[] FGDescription = {
            "terminal primary C(sp3)",
            "total secondary C(sp3)",
            "total tertiary C(sp3)",
            "total quaternary C(sp3)",
            "ring secondary C(sp3)",
            "ring tertiary C(sp3)",
            "ring quaternary C(sp3)",
            "aromatic C(sp2)",
            "unsubstituted benzene C(sp2)",
            "substituted benzene C(sp2)",
            "non-aromatic conjugated C(sp2)",
            "terminal primary C(sp2)",
            "aliphatic secondary C(sp2)",
            "aliphatic tertiary C(sp2)",
            "allenes groups",
            "terminal C(sp)",
            "non-terminal C(sp)",
            "cyanates (aliphatic)",
            "cyanates (aromatic)",
            "isocyanates (aliphatic)",
            "isocyanates (aromatic)",
            "thiocyanates (aliphatic)",
            "thiocyanates (aromatic)",
            "isothiocyanates (aliphatic)",
            "isothiocyanates (aromatic)",
            "carboxylic acids (aliphatic)",
            "carboxylic acids (aromatic)",
            "esters (aliphatic)",
            "esters (aromatic)",
            "primary amides (aliphatic)",
            "primary amides (aromatic)",
            "secondary amides (aliphatic)",
            "secondary amides (aromatic)",
            "tertiary amides (aliphatic)",
            "tertiary amides (aromatic)",
            "(thio-) carbamates (aliphatic)",
            "(thio-) carbamates (aromatic)",
            "acyl halogenides (aliphatic)",
            "acyl halogenides (aromatic)",
            "thioacids (aliphatic)",
            "thioacids (aromatic)",
            "dithioacids (aliphatic)",
            "dithioacids (aromatic)",
            "thioesters (aliphatic)",
            "thioesters (aromatic)",
            "dithioesters (aliphatic)",
            "dithioesters (aromatic)",
            "aldehydes (aliphatic)",
            "aldehydes (aromatic)",
            "ketones (aliphatic)",
            "ketones (aromatic)",
            "urea (-thio) derivatives",
            "carbonate (-thio) derivatives",
            "amidine derivatives",
            "guanidine derivatives",
            "imines (aliphatic)",
            "imines (aromatic)",
            "oximes (aliphatic)",
            "oximes (aromatic)",
            "primary amines (aliphatic)",
            "primary amines (aromatic)",
            "secondary amines (aliphatic)",
            "secondary amines (aromatic)",
            "tertiary amines (aliphatic)",
            "tertiary amines (aromatic)",
            "N hydrazines",
            "N azo-derivatives",
            "nitriles (aliphatic)",
            "nitriles (aromatic)",
            "positively charged N",
            "quaternary N",
            "hydroxylamines (aliphatic)",
            "hydroxylamines (aromatic)",
            "N-nitroso groups (aliphatic)",
            "N-nitroso groups (aromatic)",
            "nitroso groups (aliphatic)",
            "nitroso groups (aromatic)",
            "nitro groups (aliphatic)",
            "nitro groups (aromatic)",
            "imides (-thio)",
            "hydrazones",
            "hydroxyl groups",
            "aromatic hydroxyls",
            "primary alcohols",
            "secondary alcohols",
            "tertiary alcohols",
            "ethers (aliphatic)",
            "ethers (aromatic)",
            "hypohalogenides (aliphatic)",
            "hypohalogenides (aromatic)",
            "anhydrides (-thio)",
            "water molecules",
            "thiols",
            "thioketones",
            "sulfides",
            "disulfides",
            "sulfoxides",
            "sulfones",
            "sulfenic (thio-) acids",
            "sulfinic (thio-/dithio-) acids",
            "sulfonic (thio-/dithio-) acids",
            "sulfuric (thio-/dithio-) acids",
            "sulfites (thio-/dithio-)",
            "sulfonates (thio-/dithio-)",
            "sulfates (thio-/dithio-)",
            "sulfonamides (thio-/dithio-)",
            "phosphites/thiophosphites",
            "phosphates/thiophosphates",
            "phosphanes",
            "phosphonates (thio-)",
            "phosphoranes (thio-)",
            "CH2RX",
            "CHR2X",
            "CR3X",
            "R=CHX",
            "R=CRX",
            "R#CX",
            "CHRX2",
            "CR2X2",
            "R=CX2",
            "CRX3",
            "X on aromatic ring",
            "X on ring C(sp3)",
            "X on ring C(sp2)",
            "X on exo-conjugated C",
            "Aziridines",
            "Oxiranes",
            "Thiranes",
            "Azetidines",
            "Oxetanes",
            "Thioethanes",
            "Beta-Lactams",
            "Pyrrolidines",
            "Oxolanes",
            "tetrahydro-thiophenes",
            "Pyrroles",
            "Pyrazoles",
            "Imidazoles",
            "Furanes",
            "Thiophenes",
            "Oxazoles",
            "Isoxazoles",
            "Thiazoles",
            "Isothiazoles",
            "Triazoles",
            "Pyridines",
            "Pyridazines",
            "Pyrimidines",
            "Pyrazines",
            "1-3-5-Triazines",
            "1-2-4-Triazines",
            "donor atoms for H-bonds (N and O)",
            "acceptor atoms for H-bonds (N,O,F)",
            "intramolecular H-bonds (with N,O,F)"
    };


    private static String[] FGSmarts = {
            "[$([C;D1]-[#6]),$([C;D1]-[*;!#6;D1]),$([C;D2](-[#6])[*;!#6;D1]),$([C;D2](-[*;!#6;D1])[*;!#6;D1]),$([C;D3](-[#6])([*;!#6;D1])[*;!#6;D1]),$([C;D3](-[*;!#6;D1])([*;!#6;D1])[*;!#6;D1]),$([C;D4](-[#6])([*;!#6;D1])([*;!#6;D1])[*;!#6;D1]),$([C;D4](-[*;!#6;D1])([*;!#6;D1])([*;!#6;D1])[*;!#6;D1])]",  // 1
            "[$([C;D2]([#6])[#6]),$([C;D3]([#6])([#6])[*;!#6]),$([C;D4]([#6])([#6])([*;!#6])[*;!#6])]",  // 2
            "[$([C;D3]([#6])([#6])[#6]),$([C;D4]([#6])([#6])([#6])[*;!#6])]",  // 3
            "[C;D4]([#6])([#6])([#6])[#6]",  // 4
            "[$([C;D2;H2](@[#6])@[#6]),$([C;D3;H](@[#6])(@[#6])[*;!#6]),$([C;D3;H](@[#6])(@[!#6])[#6]),$([C;D4](@[#6])(@[#6])([*;!#6])[*;!#6]),$([C;D4](@[#6])(@[!#6])([#6])[*;!#6]),$([C;D4](@[!#6])(@[!#6])([#6])[#6])]",  // 5
            "[$([C;D3;H](@[#6])(@[#6])[#6]),$([C;D4](@[!#6])(@[#6])([#6])[#6]),$([C;D4](@[#6])(@[#6])([#6])[*;!#6])]",  // 6
            "[C;D4](@[#6])(@[#6])([#6])[#6]",  // 7
            "[c]",  // 8
            "[$([c;D2]1ccccc1)]",  // 9
            "[$([c;D3]1ccccc1)]",  // 10
            "[$(C=CC=*),$(C(=*)C=*),$(C(=*)[a]),$(C=[*][a])]",  // 11
            "[$([C;D1]=C),$([C;D2]([!#6;D1])=C),$([C;D3]([!#6;D1])([!#6;D1])=C)]",  // 12
            "[$([C;D2]([#6])=[#6]),$([C;D3]([!#6])([#6])=[#6])]",  // 13
            "[$([C;D3]([#6])([#6])=C)]",  // 14
            "[$([C;D2](=[#6])=[#6])]",  // 15
            "[$([C;D1]#[#6]),$([C;D2](#[#6])[*;!#6;D1])]",  // 16
            "[C;D2](#C)[C,$([*;!D1])]",  // 17
            "[$(C(#N)-O-[A])]",  // 18
            "[$(C(#N)-O-[a])]",  // 19
            "[$(C(=O)=N[A])]",  // 20
            "[$(C(=O)=N[a])]",  // 21
            "[$(C(#N)-S-[A])]",  // 22
            "[$(C(#N)-S-[a])]",  // 23
            "[$(C(=S)=N-[A])]",  // 24
            "[$(C(=S)=N-[a])]",  // 25
            "C(=O)([O;D1;!-])[A]",  // 26
            "C(=O)([O;D1;!-])[a]",  // 27
            "O=[$([C;D2]),$([C;D3]C)][O;D2][C,c]",  // 28
            "O=C([a])[O;D2][C,c]",  // 29
            "[$([C;D3](=O)([N;D1;!+])[C;A]),$([C;D2](=O)[N;D1;!+])]",  // 30
            "[$([C;D3](=O)([N;D1;!+])[a])]",  // 31
            "[$([C;D2]),$([C;D3]C)](=O)[N;D2;!+][C,c;!$(C=O)]",  // 32
            "C([a])(=O)[N;D2;!+][C,c;!$(C=O)]",  // 33
            "O=[$([C;D2]),$([C;D3]C)][N;D3;!+]([C,c;!$(C=O)])[C,c;!$(C=O)]",  // 34
            "O=C([a])[N;D3;!+]([C,c;!$(C=O)])[C,c;!$(C=O)]",  // 35
            "C(=[O,S])([$([O;D1]),$([S;D1]),$([O;D2](C)[A]),$([S;D2](C)[A])])[$([N;D1]),$([N;D2](C)[A]),$([N;D3](C)([A])[A])]",  // 36
            "[$(C(=[O,S])([O,S][a])N),$(C(=[O,S])([O,S])N[a])]",  // 37
            "[$(C(=O)([Cl,Br,F,I])[a])]",  // 38
            "[$(C(=O)([Cl,Br,F,I])[a])]",  // 39
            "[$(C(=O)([S;D1])C),$(C(=S)([O;D1;!-])C)]",  // 40
            "[$(C(=O)([S;D1])a),$(C(=S)([O;D1;!-])a)]",  // 41
            "[$(C(=S)([SH])[A])]",  // 42
            "[$(C(=S)([SH])[A])]",  // 43
            "[$(C(=S)([O;D2][C,c])C),$([C;D2](=S)([O;D2][C,c])),$(C(=O)([S;D2][C,c])C),$([C;D2](=O)([S;D2][C,c]))]",  // 44
            "[$(C(=S)([O;D2][C,c])[a]),$(C(=O)([S;D2][C,c])[a])]",  // 45
            "[$(C(=S)([S;D2][C,c])C),$([C;D2](=S)([S;D2][C,c]))]",  // 46
            "[$(C(=S)([S;D2][C,c])[a])]",  // 47
            "[$([C;D2](=O)C)]",  // 48
            "[$([C;D2](=O)[a])]",  // 49
            "[C;D3](=O)([C])[C]",  // 50
            "[$([C;D3](=O)([a])[C,c])]",  // 51
            "C(=[O,S])([$([#7;D3;!+](*)(*)*),$([#7;D2;!+](*)*),$([#7;D1;!+])])[$([#7;D3;!+](*)(*)*),$([#7;D2;!+](*)*),$([#7;D1;!+])]",  // 52
            "C(=[O,S])([$([O,S]);D2])[$([O,S]);D2]",  // 53
            "[$([C;D2]),$([C;D3][C,c])](=[$([#7;D1]),$([#7;D2]);!+])[$([#7;D1]),$([#7;D2]),$([#7;D3]);!$([#7]=*);!+]",  // 54
            "C([$([N;D1]),$([N;D2]),$([N;D3]);!$(N=*);!+])([$([N;D1]),$([N;D2]),$([N;D3]);!$(N=*);!+])=[$([N;D1]),$([N;D2]);!+]",  // 55
            "[$([N;D1;!+]),$([N;!+]C),$([N;!+][N;D2]=*)]=[$([C;D1]),$([C;D2]C),$([C;D3](C)C)]",  // 56
            "[$([N;D1;!+]),$([N;!+]C),$([N;!+][N;D2]=*)]=[$([C;D2][a]),$([C;D3]([a])[C,c]),$([C;D3]([a])N=*)]",  // 57
            "[$([C;D1]),$([C;D2]C),$([C;D3](C)C)]=[N;D2]O",  // 58
            "[$([C;D2][a]),$([C;D3]([C,c])[a])]=[N;D2]O",  // 59
            "[N;D1][$([C;A]);!$(C=[O,S])]",  // 60
            "[N;D1][a]",  // 61
            "[N;D2]([$([C;A]);!$(C=[O,S])])[$([C;A]);!$(C=[O,S])]",  // 62
            "[N;D2]([a])[$([c,C]);!$(C=[O,S])]",  // 63
            "[N;D3]([$([C;A]);!$(C=[O,S])])([$([C;A]);!$(C=[O,S])])[$([C;A]);!$(C=[O,S])]",  // 64
            "[N;D3]([a])([$([C,c]);!$(C=[O,S])])[$([C,c]);!$(C=[O,S])]",  // 65
            "[$([N;D1]),$([N;D2][C,c]),$([N;D3]([C,c])[C,c])]-[$([N;D1]),$([N;D2][C,c]),$([N;D3]([C,c])[C,c])]",  // 66
            "[$([N;D1]),$([N;D2][C,c])]=[$([N;D1]),$([N;D2][C,c])]",  // 67
            "N#CC",  // 68
            "N#C[a]",  // 69
            "[$([#7+,#7++,#7+++]);!$([N+]([O-])=O)]",  // 70
            "[N;D4]",  // 71
            "[$([N;D1;!+]),$([N;D2;!+]([!a])[!a]),$([N;D3;!+]([!a])([!a])[!a])]O",  // 72
            "[$([N;D2;!+]([*])[*]),$([N;D3;!+]([*])([*])[*])]([a])O",  // 73
            "O=[N;!$([N+])]N([C,c])[C,c]",  // 74 corrected from "O=NN([C,c])[C,c]"
            "O=NN([A,a])a",  // 75
            "[C][N;D2;!+;!++;!-;!--]=O",  // 76
            "[a][N;D2;!+;!++;!-;!--]=O",  // 77
            "[C][N+]([O-])=O",  // 78
            "[a][N+]([O-])=O",  // 79
            "[$([N;D2;!+]),$([N;D3](C)(C)[C,c])]([C;D3]=[O,S])[C;D3]=[O,S]",  // 80
            "N(=[$([C;D1]),$([C;D2][C,c]),$([C;D3]([C,c])[C,c])])[$([N;D1]),$([N;D2][C,c]),$([N;D3]([C,c])[C,c])]",  // 81
            "[O;D1;!-]A",  // 82
            "[O;D1;!-]a",  // 83
            "[O;D1;!-][C;D2;H2][C,c]",  // 84
            "[O;D1;!-][C;D3;H1]([C,c])[C,c]",  // 85
            "[O;D1;!-][C;D4]([C,c])([C,c])[C,c]",  // 86
            "[C;!$(C=[O,S]);!$(C#N)]O[C;!$(C=[O,S]);!$(C#N)]",  // 87
            "[#6;!$(C=O);!$(C#N)]O[a]",  // 88
            "[F,Cl,Br,I]O[A;!$(C=O);!$(C#N)]",  // 89
            "[F,Cl,Br,I]O[a]",  // 90
            "[C;D3](=[O,S])O[C;D3]=[O,S]",  // 91
            "[O;H2]",  // 92
            "[S;D1][C;!$(C=*);!$(C#*)]",  // 93
            "[S;D1]=C(C)C",  // 94
            "[#6;!$(C=O);!$(C=S);!$(C#*)][S;D2][#6;!$(C=O);!$(C=S);!$(C#*)]",  // 95
            "[#6;!$(C=O);!$(C=S);!$(C#*)][S;D2][S;D2][#6;!$(C=O);!$(C=S);!$(C#*)]",  // 96
            "[$([S;D3](=O)([C,c])[C,c]),$([S;D2](=O)C),$([S;D2](=O)S)]",  // 97
            "[$([S;D4](=O)(=O)([*;!a])[*;!a]),$([S;D3](=O)(=O)=C),$([S;D3](=O)(=O)=S)]",  // 98
            "[S;D2][O,S;D1]",  // 99
            "[S;D3](=[O,S])[O,S;D1]",  // 100
            "[S;D4](=[S,O])(=[S,O])([S,O;D1])[*;!S;!O]",  // 101
            "[S,O;D1]S([S,O])(=[S,O])=[S,O]",  // 102
            "[$([S;D3](=[S,O])[S,O;D2]),$([S;D3](=[S,O])([S,O;D2])[S,O;D2])]",  // 103
            "[S;D4](=[S,O])(=[S,O])([S,O;D2])[*;!S;!O]",  // 104
            "[S;D4](=[S,O])(=[S,O])([S,O;D2])[S,O;D2]",  // 105
            "[$([S;D2]),$([S;D3]=[O,S]),$([S;D4](=[O,S])=O)]([C,c])[$([N;D1;!+]),$([N;D2;!+]([*])[*]),$([N;D3;!+]([*])([*])[*])]",  // 106
            "[P;D3]([$([O,S][A,a])])([$([O,S][A,a])])[$([O,S])]",  // 107
            "[P;D4](=[O,S])([$([O,S][A,a])])([$([O,S][A,a])])[$([O,S])]",  // 108
            "[$([P;D3]([C,Cl,Br,F,I])([C,Cl,Br,F,I])[C,Cl,Br,F,I]),$([P;D2]([C,Cl,Br,F,I])[C,Cl,Br,F,I]),$([P;D1][C,Cl,Br,F,I])]",  // 109
            "[$([P;D3]),$([P;D4][C,c,F,Cl,Br,I])]([O,S])(=[O,S])[O,S]",  // 110
            "[$([#6,F,Cl,Br,I]P([#6,F,Cl,Br,I])([#6,F,Cl,Br,I])=[O,S]),$([#6,F,Cl,Br,I][P;D3]([#6,F,Cl,Br,I])=[O,S]),$([#6,F,Cl,Br,I][P;D2]=[O,S]),$([P;D1]=[O,S])]",  // 111
            "[C,c][C;D2;!R][Cl,Br,F,I]",  // 112
            "[C,c][C;D3;!R;!R2;!$(C(@[*])@[*])]([C,c])[Cl,Br,F,I]",  // 113
            "[C;D4;!R;!R2;!$(C(@[*])@[*])]([C,c])([C,c])([C,c])[Cl,Br,F,I]",  // 114
            "C=[C;D2;!R][Cl,Br,F,I]",  // 115
            "C=[C;!R;!R2;!$(C(@[*])@[*])]([C,c])[Cl,Br,F,I]",  // 116
            "C#C[Cl,Br,F,I]",  // 117
            "[C;D3;!R]([C,c])([Cl,Br,F,I])[Cl,Br,F,I]",  // 118
            "[C;D4;!R;!R2;!$(C(@[*])@[*])]([C,c])([C,c])([Cl,Br,F,I])[Cl,Br,F,I]",  // 119
            "C=[C;!R]([Cl,Br,F,I])[Cl,Br,F,I]",  // 120
            "[C;D4;!R]([C,c])([Cl,Br,F,I])([Cl,Br,F,I])[Cl,Br,F,I]",  // 121
            "a[Cl,Br,F,I]",  // 122
            "[Cl,Br,F,I][$(C(@[*])@[*]);!$(C=*)]",  // 123
            "[Cl,Br,F,I][$(C(@[*])@[*]);$(C=*)]",  // 124
            "[Cl,Br,F,I][$(C(=[*])C=[*]),$(C=[*]C=[*]),$(C=[*][a])]",  // 125
            "C1CN1",  // 126
            "C1CO1",  // 127
            "C1CS1",  // 128
            "C1C[N;!$(N@C=O)]C1",  // 129
            "C1COC1",  // 130
            "C1CSC1",  // 131
            "O=[$(C1CC2@*@*@*N12),$(C1CC2@*@*@*@*N12)]",  // 132
            "C1CC[N;!$(N@C=O)]C1",  // 133
            "C1CC[O;!$(O@C=O)]C1",  // 134
            "C1CCSC1",  // 135
            "n1cccc1",  // 136
            "n1cccn1",  // 137
            "n1ccnc1",  // 138
            "c1ccoc1",  // 139
            "c1ccsc1",  // 140
            "c1cocn1",  // 141
            "c1cnoc1",  // 142
            "c1cscn1",  // 143
            "c1cnsc1",  // 144
            "[$(c1ncnn1),$(c1cnnn1)]",  // 145
            "c1cc[$([nX3]),$([nX2])]cc1",  // 146
            "c1ccnnc1",  // 147
            "c1cncnc1",  // 148
            "c1cnccn1",  // 149
            "c1ncncn1",  // 150
            "c1cnncn1",  // 151
            "",    //"[!H0;#7,#8]",  // 152
            "",     //"[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]",  // 153
    };


    /**
     * Constructor.
     */
    public FG() {
        super();
        this.Name = FG.BlockName;
        int FGNumber = FGSmarts.length;
//        Queries = new QueryAtomContainer[FGNumber];
        Queries = new Pattern[FGNumber];
        for (int i = 0; i < FGNumber; i++) {
            try {
                Queries[i] = SmartsPattern.create(FGSmarts[i]).setPrepare(false);
            } catch (Exception e) {
                Queries[i] = null;
            }
        }
    }


    @Override
    protected final void GenerateDescriptors() {
        DescList.clear();
        int FGNumber = FGSmarts.length;
        for (int i=0; i<FGNumber; i++) {
            Add(FGNames[i], FGDescription[i]);
        }
        SetAllValues(Descriptor.MISSING_VALUE);
    }


    /**
     * Calculate descriptors for the given molecule.
     *
     * @param mol molecule to be calculated
     */
    @Override
    public void Calculate(InsilicoMolecule mol) {

        GenerateDescriptors();

        try {

            // Matcher tool
//            CustomQueryMatcher Matcher;
//            try {
//                Matcher = new CustomQueryMatcher(mol);
//            } catch (Exception e) {
//                throw new Exception("Unable to init SMARTS matcher");
//            }

            // Performs SMARTS matching
            for (int i=0; i<this.GetSize(); i++) {

                // last two groups skipped (H-donor and H-acceptor)
                if (i>=151)
                    continue;

                // Check if query has been correctly initialized
                if (Queries[i] == null)
                    throw new Exception("Unable to init SMARTS query no. " + i);

                int nmatch = 0;
//                List mappings;
                List<Mappings> mappings = new ArrayList<>();
                boolean status;
                boolean err = false;

                try {
                    status = Queries[i].matches(mol.GetStructure());
                    if (status) {
                        mappings.add(Queries[i].matchAll(mol.GetStructure()));
                        nmatch = Queries[i].matchAll(mol.GetStructure()).countUnique();
                    }
                } catch (Exception e) {
                    err = true;
                }

                // manual fix for triazoles
                if (i == 144) {
                    nmatch = nmatch / 2;
                }

                // Sets group
                if (!err)
                    this.SetByIndex(i, nmatch);
                else
                    this.SetByIndex(i, -999);

            }


            // Calculates H donor and acceptors
            IAtomContainer m = mol.GetStructure();
            int HAcc=0, HDon=0;
            int nSK = m.getAtomCount();
            for (int i=0; i<nSK; i++) {

                IAtom at = m.getAtom(i);
                int nH = 0;
                try {
                    nH = at.getImplicitHydrogenCount();
                } catch (Exception e) { }

                // H Donors: number of h linked to any N and O
                if ( (at.getSymbol().equalsIgnoreCase("O")) ||
                        (at.getSymbol().equalsIgnoreCase("N")) ) {
                    HDon += nH;
                }

                // H Acceptors
                if (at.getSymbol().equalsIgnoreCase("F"))
                    HAcc++;
                if (at.getSymbol().equalsIgnoreCase("O")) {
                    if ((6 - MoleculeUtilities.GetTotalBondOrder(at, m) - at.getFormalCharge()) > 0)
                        HAcc++;
                }
                if (at.getSymbol().equalsIgnoreCase("N")) {
                    if (MoleculeUtilities.GetTotalBondOrder(at, m) < 4)
                        HAcc++;
                    // TODO! Exclude pyrrole-like N
                }

            }
            this.SetByIndex(151, HDon);
            this.SetByIndex(152, HAcc);

        } catch (Exception e) {
            this.SetAllValues(-999);
        }

    }


    /**
     * Clones the actual descriptor block
     * @return a cloned copy of the actual object
     * @throws CloneNotSupportedException
     */
    public DescriptorBlock CreateClone() throws CloneNotSupportedException {
        FG block = new FG();
        block.CloneDetailsFrom(this);
        block.Queries = this.Queries.clone();
        return block;
    }

}
