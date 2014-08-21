/* Annotators from GATK-lite 2.3 tree, ported to current GATK base library.
*/

package bcbio.gatk.tools.walkers.annotator;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.BaseUtils;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * The GC content (# GC bases / # all bases) of the reference within 50 bp +/- this site
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
public class GCContent extends InfoFieldAnnotation implements ExperimentalAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        double content = computeGCContent(ref);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", content));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("GC"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("GC", 1, VCFHeaderLineType.Integer, "GC content within 20 bp +/- the variant")); }

    public boolean useZeroQualityReads() { return false; }

    private static double computeGCContent(ReferenceContext ref) {
        int gc = 0, at = 0;

        for ( byte base : ref.getBases() ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            if ( baseIndex == 1 || baseIndex == 2 )
                gc++;
            else if ( baseIndex == 0 || baseIndex == 3 )
                at++;
            else
                ; // ignore
        }

        int sum = gc + at;
        return (100.0*gc) / (sum == 0 ? 1 : sum);
     }
}
