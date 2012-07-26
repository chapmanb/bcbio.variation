/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package bcbio.variation.annotate;

import cern.jet.math.Arithmetic;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


/**
 * Calculate Most Probable Genotype similar to NHGRI's bam2mpg (research.nhgri.nih.gov/software/bam2mpg/overview.shtml)
 *   This also calculates the raw likelihoods for each genotype so that they can be combined across datasets
 * Note that this is currently not calculated for indels or for multi-allelic sites.
 * This code was adapted from FisherStrand by Justin Zook at NIST (jzook at nist.gov)
 */
public class MostProbableGenotype extends InfoFieldAnnotation implements ExperimentalAnnotation {
    private static final String MPG = "MPG";

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
//        if ( ! vc.isVariant() || vc.isFiltered() )
//            return null;

        double[] MPGs;

        if (vc.isBiallelic() && vc.isSNP())
            MPGs = calcSNPMPGs(stratifiedContexts, vc.getReference(), vc.getAlternateAllele(0));
/*        else if (vc.isIndel() || vc.isMixed()) {
            table = getIndelContingencyTable(stratifiedContexts, vc);
            if (table == null)
                return null;
        }
*/        else
            return null;

        if ( MPGs == null )
            return null;

        // output MPG ratios and raw likelihoods
        Map<String, Object> map = new HashMap<String, Object>();
            map.put("MPGHomRef",MPGs[0]);
            map.put("MPGHetRef",MPGs[1]);
            map.put("MPGHomVar",MPGs[2]);
            map.put("MPGLikHomRef",MPGs[3]);
            map.put("MPGLikHetRef",MPGs[4]);
            map.put("MPGLikHomVar",MPGs[5]);
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList("MPGHomRef","MPGHetRef","MPGHomVar","MPGLikHomRef","MPGLikHetRef","MPGLikHomVar"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MPGHomRef", 1, VCFHeaderLineType.Float, "Most Probable Genotype likelihood ratio for Hom Ref"),
            new VCFInfoHeaderLine("MPGHetRef", 1, VCFHeaderLineType.Float, "Most Probable Genotype likelihood ratio for Het Ref"),
            new VCFInfoHeaderLine("MPGHomVar", 1, VCFHeaderLineType.Float, "Most Probable Genotype likelihood ratio for Hom Var"),
            new VCFInfoHeaderLine("MPGLikHomRef", 1, VCFHeaderLineType.Float, "Raw likelihood for Hom Ref"),
            new VCFInfoHeaderLine("MPGLikHetRef", 1, VCFHeaderLineType.Float, "Raw likelihood for Het Ref"),
            new VCFInfoHeaderLine("MPGLikHomVar", 1, VCFHeaderLineType.Float, "Raw likelihood for Hom Var")); }



    /**
     Calculate confidence intervals for both strands
     */
    private static double[] calcSNPMPGs(Map<String, AlignmentContext> stratifiedContexts, Allele ref, Allele alt) {
        int steps = 2; // number of divisions between 0 and 1 to use for allele balance when calculating probabilities
        double[] pBoth = new double[steps+1]; // log sum of probabilities for both strands
	double[] MPGs = new double[6]; 
	Arrays.fill(pBoth, 0);
	Arrays.fill(MPGs, 0);

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for (PileupElement p : sample.getValue().getBasePileup()) {
                if ( p.isDeletion()  ) // ignore deletions 
                    continue;

                //if ( p.getRead().getMappingQuality() < 20 || p.getQual() < 10 )
                //    continue; // todo -- fixme, should take filtered context!
 

                Allele base = Allele.create(p.getBase(), false);
                 boolean matchesRef = ref.equals(base, true);
                boolean matchesAlt = alt.equals(base, true);
		double eps = Math.pow(10.0,-1.0*(p.getQual()/10.0)); // convert base quality score to error rate
		int vAB=0;
	        double[] p1 = new double[steps+1]; // log probability for this base (or log(p(b|eps,AB))), similar to ContEst with f=0

		//Arrays.fill(p1, 0);
                if ( matchesRef ) {
			vAB=0;
			while (vAB <= steps) {
				p1[vAB] = Math.log10( ((double) vAB/(double) steps)*(1-eps) + (1-((double) vAB/(double) steps))*eps );
				vAB++;
			}
		} else { // if ( matchesAlt ) {
			vAB=0;
			while (vAB <= steps) {
				p1[vAB] = Math.log10( 1-(((double) vAB/(double) steps)*(1-eps) + (1-((double) vAB/(double) steps))*eps) );
				vAB++;
			}
		} 
		// add log probabilities for both strands for each step
			vAB=0;
			while (vAB <= steps) {
				pBoth[vAB] += p1[vAB];
				vAB++;
			}
 		// add log probabilities for negative strands for each step
/*               if ( p.getRead().getReadNegativeStrandFlag() ) {
			vAB=0;
			while (vAB <= steps) {
				pRev[vAB] += p1[vAB];
				vAB++;
			}
		} else {  	// add log probabilities for positive strands for each step
			vAB=0;
			while (vAB <= steps) {
				pFor[vAB] += p1[vAB];
				vAB++;
			}
		}
*/

            }
        }
	int vAB=0;
//			double listSum = 0;
/*	while (vAB <= steps) {
		if (pBoth[vAB] != 0) pBoth[vAB] = Math.pow(10,pBoth[vAB]);
//		if (pRev[vAB] != 0) pRev[vAB] = Math.pow(10,pRev[vAB]);
//		if (pFor[vAB] != 0) pFor[vAB] = Math.pow(10,pFor[vAB]);
		vAB++;
	}
*/
	// Calc genotype likelihood ratios for HomRef, HetRef, and HomVar
	if (pBoth[0]<pBoth[1]) { //HomRef
	    MPGs[0] = (pBoth[2] - pBoth[1]);
	} else {
	    MPGs[0] = (pBoth[2] - pBoth[0]);
	}
	if (pBoth[0]<pBoth[2]) { //HetRef
	    MPGs[1] = (pBoth[1] - pBoth[2]);
	} else {
	    MPGs[1] = (pBoth[1] - pBoth[0]);
	}
	if (pBoth[2]<pBoth[1]) { //HomVar
	    MPGs[2] = (pBoth[0] - pBoth[1]);
	} else {
	    MPGs[2] = (pBoth[0] - pBoth[2]);
	}

	// Also return raw likelihoods for HomRef, HetRef, and HomVar so that they can be combined across datasets
	MPGs[3] = -1*(pBoth[2]);
	MPGs[4] = -1*(pBoth[1]);
	MPGs[5] = -1*(pBoth[0]);

        return MPGs;
    }



    /**
	TODO: enable MPG calculation for indels
     */
    private static int[][] getIndelContingencyTable(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        final double INDEL_LIKELIHOOD_THRESH = 0.3;
        final HashMap<PileupElement,LinkedHashMap<Allele,Double>> indelLikelihoodMap = IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap();

        if (indelLikelihoodMap == null)
            return null;
        
        int[][] table = new int[2][2];

        for ( String sample : stratifiedContexts.keySet() ) {
            final AlignmentContext context = stratifiedContexts.get(sample);
            if ( context == null ) 
                continue;

            ReadBackedPileup pileup = context.getBasePileup();
            if (pileup == null)
                 continue;

            for (final PileupElement p: pileup) {
                if ( p.getRead().getMappingQuality() < 20)
                    continue;
                if (indelLikelihoodMap.containsKey(p)) {
                    // to classify a pileup element as ref or alt, we look at the likelihood associated with the allele associated to this element.
                    // A pileup element then has a list of pairs of form (Allele, likelihood of this allele).
                    // To classify a pileup element as Ref or Alt, we look at the likelihood of corresponding alleles.
                    // If likelihood of ref allele > highest likelihood of all alt alleles  + epsilon, then this pileup element is "ref"
                    // otherwise  if highest alt allele likelihood is > ref likelihood + epsilon, then this pileup element it "alt"
                    // retrieve likelihood information corresponding to this read
                    LinkedHashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                    // by design, first element in LinkedHashMap was ref allele
                    boolean isFW = !p.getRead().getReadNegativeStrandFlag();

                    double refLikelihood=0.0, altLikelihood=Double.NEGATIVE_INFINITY;

                    for (Allele a : el.keySet()) {

                        if (a.isReference())
                            refLikelihood =el.get(a);
                        else {
                            double like = el.get(a);
                            if (like >= altLikelihood)
                                altLikelihood = like;
                        }
                    }

                    boolean matchesRef = (refLikelihood > (altLikelihood + INDEL_LIKELIHOOD_THRESH));
                    boolean matchesAlt = (altLikelihood > (refLikelihood + INDEL_LIKELIHOOD_THRESH));
                    if ( matchesRef || matchesAlt ) {
                        int row = matchesRef ? 0 : 1;
                        int column = isFW ? 0 : 1;

                         table[row][column]++;
                    }


                }
            }
        }

        return table;
    }
}
