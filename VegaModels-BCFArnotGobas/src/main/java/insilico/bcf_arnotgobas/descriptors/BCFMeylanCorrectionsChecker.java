package insilico.bcf_arnotgobas.descriptors;

import insilico.core.alerts.Alert;
import insilico.core.alerts.AlertBlockFromSMARTS;
import insilico.core.alerts.AlertEncoding;
import insilico.core.alerts.AlertList;
import insilico.core.alerts.iAlertBlock;
import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;

/**
 *
 * @author amanganaro
 */
public class BCFMeylanCorrectionsChecker extends AlertBlockFromSMARTS implements iAlertBlock  {

    private double[] CorrFactor;
    private String[] SMARTSFragments;
    private Pattern[] SA;
    
    private double Correction;
    private double CurLogP;
    private boolean TinAndMercury;
    
    
    public BCFMeylanCorrectionsChecker() throws Exception {
        super(-1, "Meylan BCF Correction Fragments"); // id = -1 as it is not needed (it is not a public alert block)
    }
    
    
    @Override
    protected void BuildSAList() throws InitFailureException {

        CorrFactor = new double[12];
        SMARTSFragments = new String[12];
        
        Alert curSA;
        int i=0;
        
        SMARTSFragments[i] = "[c]C(=O)[C,c]";
        CorrFactor[i] = -0.5851;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Ketone (aromatic connection)");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "P(=O)([O;D2])([O;D2])(O)";
        CorrFactor[i] = -0.8250;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Phosphate ester");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = ""; // done by code
        CorrFactor[i] = 0.5860;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Multi-halogenated biphenyl");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = ""; // done by code
        CorrFactor[i] = 0.5860;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Multi-halogenated PAH");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = "[a][C;D2,D3][OH]";
        CorrFactor[i] = -0.2556;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Aromatic ring-CH-OH");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = "c1ncncn1";
        CorrFactor[i] = -0.5169;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Sym-triazine");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "c1cc([O;D1])c(C([C;D1])([C;D1])[C;D1])cc1";
        CorrFactor[i] = -0.2220;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Tert-butyl ortho-phenol");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "[c;R1]1[c;R1][c;R1]c2c([c;R1]1)[c;R1][c;R1]c3[c;R1][c;R1][c;R1][c;R1]c23";
        CorrFactor[i] = 0.6609;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Phenanthrene");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = "C1CC1C(=O)O[A]";
        CorrFactor[i] = -1.2591;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Cyclopropyl-C(=O)-O-");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "[C;D2;!R]";
        CorrFactor[i] = 0; // factor depends on LogP value
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Alkyl chain (8+ carbons)");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;
        
        SMARTSFragments[i] = "[S;D2;!R]-[S;D2;!R]";
        CorrFactor[i] = -1.3404;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Disulfide");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);
        i++;

        SMARTSFragments[i] = "[$([Sn,Hg]);!$([Sn,Hg][N,n])]";
        CorrFactor[i] = 1.4;
        curSA = new Alert(BlockIndex, AlertEncoding.BuildAlertId(BlockIndex, (i+1)));
        curSA.setName("MeyBCFCor" + (i+1));
        curSA.setDescription("Meylan BCF correction fragment n. " + (i+1) + ": " +
                "Tin or Mercury compound");
        curSA.setNumericProperty("COEF", CorrFactor[i]);
        Alerts.add(curSA);

    }    
    
    
    @Override
    protected void InitSMARTS() throws InitFailureException {

        try {

//            SA = new QueryAtomContainer[SMARTSFragments.length];

            SA = new Pattern[SMARTSFragments.length];

            for (int i=0; i<SMARTSFragments.length; i++)  {
                if (SMARTSFragments[i].isEmpty())
                    SA[i] = null;
                else
                    SA[i] = SmartsPattern.create(SMARTSFragments[i], DefaultChemObjectBuilder.getInstance());
            }

        } catch (Exception e) {
            throw new InitFailureException("Unable to initialize SMARTS");
        }
    }
    
    
    @Override
    protected AlertList CalculateSAMatches() throws GenericFailureException {
        
        AlertList Res = new AlertList();
        Correction = 0;
        TinAndMercury = false;
        
        try {

            for (int i=0; i<SA.length; i++) {
                
                if (i == 2) {
                    if (IsMultiHalogenatedPAHorBiphenyl()) {
                        Res.add((Alert)Alerts.get(i).clone());
                        Correction += CorrFactor[2];
                    }
                    continue;
                }
                
                // skips - already done in #2
                if (i == 3) 
                    continue;
                
                if (SA[i].matches(CurMol.GetStructure())) {
                    
                    // Alkyl chain

                    if (i == 9) {
                        int nMatches = SA[i].matchAll(CurMol.GetStructure()).countUnique();
                        if (nMatches >= 7) {
                            if ((CurLogP > 4.0) && (CurLogP < 7.0)) {
                                Res.add((Alert)Alerts.get(i).clone());
                                Correction += -1.3743;
                            }
                            if ((CurLogP >= 7.0) && (CurLogP <= 10.0)) {
                                Res.add((Alert)Alerts.get(i).clone());
                                Correction += -0.5965;
                            }
                        } 

                        continue;
                    }

                    if (i == 11)
                        TinAndMercury = true;

                    // Other SMARTS based corrections
                    Res.add((Alert)Alerts.get(i).clone());
                    Correction += CorrFactor[i];
                }
            }
            
        } catch (CloneNotSupportedException | InvalidMoleculeException e) {
            return null;
        }
        
        return Res; 
    }
    
    
    public double GetCorrection() {
        return Correction;
    }    
    
    
    public void SetLogP(double logP) {
        this.CurLogP = logP;
    }
    
    
    public boolean GetTinAndMercury() {
        return this.TinAndMercury;
    }

    
    private boolean IsMultiHalogenatedPAHorBiphenyl() throws InvalidMoleculeException {
        
        IRingSet singleRings = CurMol.GetSSSR();
        int nRings = singleRings.getAtomContainerCount();

        for (int i=0; i<nRings; i++) {
            IRing curRing = (IRing) singleRings.getAtomContainer(i);
            if (IsBenzene(curRing)) {
            
                boolean PAHorBiphenyl = false;
                int nHalo = 0;
                
                // Checks if PAH
                IRingSet ConnRings = singleRings.getConnectedRings(curRing);
                if (ConnRings != null) {
                    for (int k=0; k<ConnRings.getAtomContainerCount(); k++) {
                        if (IsBenzene((IRing) ConnRings.getAtomContainer(k))) {
                            PAHorBiphenyl = true; break;
                        }
                    }
                }
                
                // Checks if biphenyl and for halo-substituents
                for (int k=0; k<curRing.getAtomCount(); k++) {
                    IAtom at = curRing.getAtom(k);
                    List<IAtom> ConnAt = CurMol.GetStructure().getConnectedAtomsList(at);
                    for (int j=0; j< ConnAt.size(); j++) {
                        if (curRing.contains(ConnAt.get(j)))
                            continue;
                        if (ConnAt.get(j).getAtomicNumber() == 9) {
                            nHalo++; continue;
                        }
                        if (ConnAt.get(j).getAtomicNumber() == 17) {
                            nHalo++; continue;
                        }
                        if (ConnAt.get(j).getAtomicNumber() == 35) {
                            nHalo++; continue;
                        }
                        if (ConnAt.get(j).getAtomicNumber() == 53) {
                            nHalo++; continue;
                        }
                        if (ConnAt.get(j).getAtomicNumber() == 6) {
                            // Checks if biphenyl
                            if (singleRings.contains(ConnAt.get(j))) {
                                IRingSet buf = singleRings.getRings(ConnAt.get(j));
                                if (buf.getAtomContainerCount() != 1)
                                    continue;
                                if (IsBenzene((IRing) buf.getAtomContainer(0))) {
                                    PAHorBiphenyl = true; continue;
                                }
                            }
                        }
                    }
                }
                
                if ((PAHorBiphenyl) && (nHalo > 1))
                    return true;
                
            }
        }
        
        return false;
    }
    
    
    private boolean IsBenzene(IRing r) {
        
        boolean answ = true;
        
        for (int i=0; i<r.getAtomCount(); i++) {
            IAtom at = r.getAtom(i);
            if (at.getAtomicNumber() != 6) {
                answ = false; break;
            }
            if (at.getFlag(CDKConstants.ISAROMATIC) == false) {
                answ = false; break;
            }
        }
        
        return answ;
    }
}
