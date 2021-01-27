package insilico.bcf_arnotgobas.descriptors;

import insilico.core.ad.reasoning.DescriptorAnalysis;

/**
 *
 * @author Alberto Manganaro (a.manganaro@kode-solutions.net)
 */
public class DescriptorAnalysisBCFMeylan extends DescriptorAnalysis {

    private final static String DescriptorName = "LogP (Meylan)";
    private final static String ExperimentalName = "LogBCF";
    private final static String LongDescription = "LogP is directly correlated to the logBCF value.";
    
    
    public DescriptorAnalysisBCFMeylan(int DescriptorIdx) {
        super(DescriptorIdx);
    }
    
    @Override
    public String getDescriptorName() {
        return DescriptorName;
    }

    @Override
    public String getExpName() {
        return ExperimentalName;
    }

    @Override
    public String getDescription() {
        return LongDescription;
    }
    
}
