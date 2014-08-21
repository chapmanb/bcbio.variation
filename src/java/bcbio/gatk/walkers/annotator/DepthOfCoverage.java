/* Annotators from GATK-lite 2.3 tree, ported to current GATK base library.
*/

package bcbio.gatk.tools.walkers.annotator;

import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Total (unfiltered) depth over all samples.
 *
 * While the sample-level (FORMAT) DP field describes the total depth of reads that passed the Unified Genotyper's
 * internal quality control metrics (like MAPQ > 17, for example), the INFO field DP represents the unfiltered depth
 * over all samples.  Note though that the DP is affected by downsampling (-dcov), so the max value one can obtain for
 * N samples with -dcov D is N * D
 */
public class DepthOfCoverage extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {

        int depth = 0;
        if (stratifiedContexts != null) {
            if ( stratifiedContexts.size() == 0 )
                return null;

            for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() )
                depth += sample.getValue().getBasePileup().depthOfCoverage();
        }
        else if (perReadAlleleLikelihoodMap != null) {
            if ( perReadAlleleLikelihoodMap.size() == 0 )
                return null;

            for (PerReadAlleleLikelihoodMap maps : perReadAlleleLikelihoodMap.values() ) {
                for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                    final GATKSAMRecord read = el.getKey();
                    depth += 1;
                }
            }
        }
        else
            return null;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", depth));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.DEPTH_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
