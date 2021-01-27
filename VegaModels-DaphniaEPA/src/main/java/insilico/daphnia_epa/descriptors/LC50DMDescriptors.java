package insilico.daphnia_epa.descriptors;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author amanganaro
 */
public class LC50DMDescriptors {

    private double xc4;
    private double ssi; // TODO!!!!!!
    private int n_CdS2N;
    private int n_AN;
    private int n_NP;
    private int n_SdOdO;

    Pattern[] q;
    
    
    public LC50DMDescriptors() throws InitFailureException {
        
        q = new Pattern[4];
        try {
            q[0] = SmartsPattern.create("N-C(=S)-N", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            q[1] = SmartsPattern.create("n", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            q[2] = SmartsPattern.create("P-N(-*)-*", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
            q[3] = SmartsPattern.create("a-S(=O)(=O)-*", DefaultChemObjectBuilder.getInstance()).setPrepare(false);
        } catch (Exception ex) {
            throw new InitFailureException("unable to init SMARTS parsers");
        }
        
    }
    
    
    public void Calculate(InsilicoMolecule Mol) throws GenericFailureException {
        
        IAtomContainer mol;
        try {
            mol = Mol.GetStructure();
        } catch (InvalidMoleculeException e) {
            throw new GenericFailureException(e.getMessage());
        }
        
        int nSK = mol.getAtomCount();
        
        // Chi cluster descriptor
        double DegreeMat[] = new double[nSK];
        for (int i=0; i<nSK; i++) 
            DegreeMat[i] = mol.getConnectedBondsCount(mol.getAtom(i));
        LinkedList clusters = FindClusters(mol, 4);
        xc4 = CalcChi(clusters, DegreeMat);
        
        
        try {

            List mappings = null;
            boolean status = false;
            
//            CustomQueryMatcher Matcher;
//            Matcher = new CustomQueryMatcher(Mol);
            
            n_CdS2N = 0;
            status = q[0].matches(mol);
            if (status) {
                n_CdS2N = q[0].matchAll(mol).countUnique();
            }
            
            n_AN = 0;
            status = q[1].matches(mol);
            if (status) {
                n_AN = q[1].matchAll(mol).countUnique();
            }
            
            n_NP = 0;
            status = q[2].matches(mol);
            if (status) {
                n_NP =  q[2].matchAll(mol).countUnique();
            }
            
            n_SdOdO = 0;
            status = q[3].matches(mol);
            if (status) {
                n_SdOdO = q[3].matchAll(mol).countUnique();
            }
            
        } catch (Exception e) {
            throw new GenericFailureException("unable to perform matching - " + e.getMessage());
        }
                
        
    }
     
    
    private double CalcChi(LinkedList MasterList, double[] D) {
        // in this method Chi is calculated for an LinkedList
        // which contains one or more array lists with the atom numbers

        double CHI = 0;

        //System.out.println(MasterList.size());

        for (int i = 0; i <= MasterList.size() - 1; i++) {

                LinkedList<Integer> al = (LinkedList) MasterList.get(i);

                double product = 1;

                for (int j = 0; j <= al.size() - 1; j++) {			
                        int J = al.get(j);
                        product *= D[J];
                }
                CHI += Math.pow(product, -0.5);

        }
        return CHI;
    }
    

    private static LinkedList FindClusters(IAtomContainer m, int M) {

            //finds clusters of length 3 or 4 when have atom with 3 or 4 atoms attached


            LinkedList MasterList = new LinkedList();

            for (int I = 0; I <= m.getAtomCount() - 1; I++) { // start of atom 1-
                    // need to
                    // start process with each
                    // atom to get all possible
                    // paths

                    List ca = m.getConnectedAtomsList(m.getAtom(I));

                    if (ca.size() == M) {
                            LinkedList al = new LinkedList();

                            al.add(m.indexOf(m.getAtom(I)));
                            for (int j = 0; j <= M - 1; j++) {					
                                    al.add(m.indexOf(((IAtom) ca.get(j))));
                            }
                            MasterList.add(al);
                    } else if (ca.size() > M) { 

                            LinkedList mlist=FindAllPossibleCombos(M,ca.size());

                            for (int i=0;i<=mlist.size()-1;i++) {
                                    LinkedList al = new LinkedList();
                                    al.add(m.indexOf(m.getAtom(I)));
                                    List l= Parse((String)mlist.get(i),"\t");
                                    for (int j=0;j<=l.size()-1;j++) {
                                            String snum=(String)l.get(j);
                                            int num= Integer.parseInt(snum);
                                            al.add(m.indexOf(((IAtom) ca.get(num))));
                                    }
                                    MasterList.add(al);

                            }



                    }

            } 

            return MasterList;

    
        }
	
	
        private static LinkedList FindAllPossibleCombos(int M, int NBound) {
            String s1="",s2="",s3="",s4="";

            LinkedList mlist=new LinkedList();

            for (int i=0;i<=NBound-1;i++) {

                    s1=i+"";

                    for (int j=i+1;j<=NBound-1;j++) {

                            s2=s1+"\t"+j;

                            for (int k=j+1;k<=NBound-1;k++) {

                                    s3=s2+"\t"+k;

                                    if (M==3) {
                                            mlist.add(s3);
                                    } else if (M==4) {

                                            for (int l=k+1;l<=NBound-1;l++) {
                                                    s4=s3+"\t"+l;
                                                    mlist.add(s4);
                                            } // end l for loop

                                    }// end M==4 else if


                            } // end k for loop

                    } // end j for loop

            } // end i for loop

            return mlist;

	}
	

      private static LinkedList Parse(String Line, String Delimiter) {
        // parses a delimited string into a list

        LinkedList myList = new LinkedList();

        int tabpos = 1;

        while (tabpos > -1) {
          tabpos = Line.indexOf(Delimiter);

          if (tabpos > 0) {
            myList.add(Line.substring(0, tabpos));
            Line = Line.substring(tabpos + Delimiter.length(), Line.length());
          } else if (tabpos == 0) {
            myList.add("");
            Line = Line.substring(tabpos + Delimiter.length(), Line.length());
          } else {
            myList.add(Line.trim());
          }
        }

        return myList;
      }    

    /**
     * @return the xc4
     */
    public double getXc4() {
        return xc4;
    }

    /**
     * @return the ssi
     */
    public double getSsi() {
        return ssi;
    }

    /**
     * @return the n_CdS2N
     */
    public int getN_CdS2N() {
        return n_CdS2N;
    }

    /**
     * @return the n_AN
     */
    public int getN_AN() {
        return n_AN;
    }

    /**
     * @return the n_NP
     */
    public int getN_NP() {
        return n_NP;
    }

    /**
     * @return the n_SdOdO
     */
    public int getN_SdOdO() {
        return n_SdOdO;
    }
    
}
