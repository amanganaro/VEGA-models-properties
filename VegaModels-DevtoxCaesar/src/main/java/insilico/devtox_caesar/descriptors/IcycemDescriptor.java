package insilico.devtox_caesar.descriptors;

import insilico.core.descriptor.Descriptor;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import org.openscience.cdk.RingSet;

/**
 *
 * @author Alberto Manganaro <a.manganaro@kode-solutions.net>
 */
public class IcycemDescriptor {

    public static double Calculate(InsilicoMolecule Mol) {
        
        try {
            
            double icyce;
            double icycem;

            RingSet rings = Mol.GetAllRings();
            int C = rings.getAtomContainerCount();
            int A = Mol.GetStructure().getAtomCount();

            if (C>0) {

                int n0=0;		
                int n1=0;

                for (int i = 0; i < rings.getAtomContainerCount(); i++) {
                    int ringSizei=rings.getAtomContainer(i).getAtomCount();
                    n1+=ringSizei;
                    n0+=A-ringSizei;			
                }

                if (n0>0) 
                    icyce=A*C*Log(2,A*C)-n1*Log(2,n1)-n0*Log(2,n0);
                else 
                    icyce=A*C*Log(2,A*C)-n1*Log(2,n1);

                icycem=icyce/(A*C);

            } else {
                    
                icycem=0;
                
            }

            return icycem;
        
        } catch (InvalidMoleculeException e) {
            return Descriptor.MISSING_VALUE;
        }
    }
    
    
    private static double Log(int base,double x) {
        double Logbx = Math.log10(x)/ Math.log10((double)base);
        return Logbx;
   }
    
}
