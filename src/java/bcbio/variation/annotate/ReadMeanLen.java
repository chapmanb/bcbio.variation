package bcbio.variation.annotate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Calculates parameters related to artifacts of incorrect local alignments
 * Mean number of aligned bases for reads covering this position - lower numbers indicative of local alignment difficulty
 * written by Justin Zook - 7/5/12
 */
public class ReadMeanLen extends InfoFieldAnnotation implements ExperimentalAnnotation {

    public List<String> getKeyNames() {
        return Arrays.asList("ReadMeanLen");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("ReadMeanLen", 1, VCFHeaderLineType.Float, "Mean number of aligned bases for reads - low number indicate possible mis-alignments"));
    }

  public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, final Map<String, PerReadAlleleLikelihoodMap> stratifiedLikelihoodMap) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        double readLenSum = 0;
		double numReads = 0;
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for ( PileupElement p : sample.getValue().getBasePileup() )	{
        	        int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, 0, 0);
        	        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips(p.getRead());

			readLenSum+=((double) numAlignedBases);

			numReads+=1;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.01f", readLenSum/numReads));
        return map;
    }



    int getNumClippedBasesAtStart(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < 20)
                numStartClippedBases++;
            else
                break;

        }

        return numStartClippedBases;
    }

    int getNumAlignedBases(SAMRecord read) {
        return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    int getNumClippedBasesAtEnd(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < 20)
                numEndClippedBases++;
            else
                break;
        }


        return numEndClippedBases;
    }

    int getOffsetFromClippedReadStart(SAMRecord read, int offset) {
        return offset - getNumClippedBasesAtStart(read);
    }
}



