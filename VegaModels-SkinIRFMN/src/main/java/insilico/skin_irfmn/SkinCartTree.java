package insilico.skin_irfmn;

/**
 *
 * @author Alberto
 */
public class SkinCartTree {

    public static boolean Predict(double nDB, double IC1, double SIC2, double GATS8s,
            double nRCHO, double CATS2D_01_NL, double F04C_O, double MLOGP) {
        
        if (MLOGP < 5.038) {
            
            if (nRCHO < 0.5) {
                
                if (nDB < 3.5) {
                    
                    if (GATS8s >= 2.6) {
                        return false;
                    } else {
                        
                        if (IC1 < 3.264) {
                            
                            if (CATS2D_01_NL >= 0.5) {
                                return false;
                            } else {
                                
                                if (nDB < 1.5) {
                                    
                                    if (SIC2 < 0.842) {
                                        
                                        if (SIC2 >= 0.681) {
                                            return false;
                                        } else {
                                            
                                            if (IC1 < 1.916) {
                                                return false;
                                            } else {
                                                return true;
                                            }
                                            
                                        }
                                        
                                    } else {
                                        return true;
                                    }
                                    
                                } else {
                                    
                                    if (F04C_O >= 2.5) {
                                        
                                        if (MLOGP < 2.806) {
                                            
                                            if (SIC2 >= 0.7735) {
                                                return false;
                                            } else {
                                                
                                                if (GATS8s >= 0.6325) {
                                                    return false;
                                                } else {
                                                    return true;
                                                }
                                                
                                            }
                                            
                                        } else {
                                            return true;
                                        }
                                        
                                    } else {
                                        return true;
                                    }
                                    
                                }
                                
                            }
                            
                        } else {
                            return true;
                        }
                        
                    }
                    
                } else {
                    return true;
                }
                
            } else {
                return true;
            }
            
        } else {
            return true;
        }
        
    }
     
}
