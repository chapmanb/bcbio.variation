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
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


/**
 * Calculates 95% confidence intervals of allele balance 
 * Note that this may not be
 * calculated for certain complex indel cases or for multi-allelic sites.
 * This code was adapted from FisherStrand by Justin Zook at NIST
 */
public class AlleleBalanceConfidenceInterval extends InfoFieldAnnotation implements ExperimentalAnnotation {
    private static final String ABCI = "ABCI";

  public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, final Map<String, PerReadAlleleLikelihoodMap> stratifiedLikelihoodMap) {
//        if ( ! vc.isVariant() || vc.isFiltered() )
//            return null;

        double[] CIs;

        if (vc.isBiallelic() && vc.isSNP())
            CIs = calcSNPCIs(stratifiedContexts, vc.getReference(), vc.getAlternateAllele(0));
/*        else if (vc.isIndel() || vc.isMixed()) {
            table = getIndelContingencyTable(stratifiedContexts, vc, stratifiedLikelihoodMap);
            if (table == null)
                return null;
        }
*/        else
            return null;

        if ( CIs == null )
            return null;

        // output AB CIs as (lower 95%, 50%, upper 95%) for both strands combined
        Map<String, Object> map = new HashMap<String, Object>();
            map.put("ABCI5",CIs[0]);
            map.put("ABCI50",CIs[1]);
            map.put("ABCI95",CIs[2]);
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList("ABCI5","ABCI50","ABCI95"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ABCI5", 1, VCFHeaderLineType.Float, "Allele Balance lower 5% confidence interval"),
            new VCFInfoHeaderLine("ABCI50", 1, VCFHeaderLineType.Float, "Allele Balance middle 50% confidence interval"),
            new VCFInfoHeaderLine("ABCI95", 1, VCFHeaderLineType.Float, "Allele Balance upper 5% confidence interval")); }



    /**
     Calculate confidence intervals for both strands
     */
    private static double[] calcSNPCIs(Map<String, AlignmentContext> stratifiedContexts, Allele ref, Allele alt) {
        int steps = 100; // number of divisions between 0 and 1 to use for allele balance when calculating probabilities
        double[] pBoth = new double[steps+1]; // log sum of probabilities for both strands
	double[] CIs = new double[3]; 
	Arrays.fill(pBoth, 0);
	Arrays.fill(CIs, 0);

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for (PileupElement p : sample.getValue().getBasePileup()) {
                if ( p.isDeletion()  ) // ignore deletions 
                    continue;

                //if ( p.getRead().getMappingQuality() < 20 || p.getQual() < 10 )
                //    continue; // todo -- fixme, should take filtered context!
                //if ( pBoth[50] != 0 )
                //    continue; // troubleshoot


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
 /*              } else {
			vAB=0;
			while (vAB <= steps) {
				p1[vAB] = Math.log10( eps/3 );
				vAB++;
			} */
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
	while (vAB <= steps) {
		if (pBoth[vAB] != 0) pBoth[vAB] = Math.pow(10,pBoth[vAB]);
//		if (pRev[vAB] != 0) pRev[vAB] = Math.pow(10,pRev[vAB]);
//		if (pFor[vAB] != 0) pFor[vAB] = Math.pow(10,pFor[vAB]);
		vAB++;
	}
	// find confidence intervals
//System.out.println("before calcCIs\n");
//	calcCIsFromProb(CIs, pFor, steps, 1);
//System.out.println("after 1 calcCIs\n");
//	calcCIsFromProb(CIs, pRev, steps, 2);
	calcCIsFromProb(CIs, pBoth, steps, 1);
//Arrays.fill(CIs, 0);
        return CIs;
    }

    /**
	Calculate 5, 50, and 95% confidence levels from list of probabilities
     */
    private static void calcCIsFromProb(double[] CIs, double[] pList, int steps, int listNum) {
		double listSum = 0;
		int vAB=0;
		while (vAB <= steps) {
			listSum += pList[vAB];
/*			if (vAB == 0) CIs[listNum*3-3] = listSum;
			if (vAB == 67) CIs[listNum*3-2] = listSum;
			if (vAB == 100) CIs[listNum*3-1] = listSum;
*/			vAB++;
		}

		if (listSum <= 0) {
			CIs[listNum*3-3] = 2;
			CIs[listNum*3-2] = 2;
			CIs[listNum*3-1] = 2;
		} else {
		double tempSum = 0;
		double ci5 = listSum * 0.05;
		double ci50 = listSum * 0.5;
		double ci95 = listSum * 0.95;
		vAB=0;
		while (tempSum <= ci5) {
			tempSum += pList[vAB];
			vAB++;
		}
		if (vAB == 0) {
			CIs[listNum*3-3] = 0;
		} else {
			CIs[listNum*3-3] = (((double) vAB - 1) - (tempSum - ci5)/pList[vAB-1])/((double) steps);
		}
		
		while (tempSum <= ci50) {
			tempSum += pList[vAB];
			vAB++;
		}
		if (vAB == 0) {
			CIs[listNum*3-2] = 0;
		} else {
			CIs[listNum*3-2] = (((double) vAB - 1) - (tempSum - ci50)/pList[vAB-1])/((double) steps);
		}

		while (tempSum <= ci95) {
			tempSum += pList[vAB];
			vAB++;
		}
		if (vAB == 0) {
			CIs[listNum*3-1] = 0;
		} else {
			CIs[listNum*3-1] = (((double) vAB - 1) - (tempSum - ci95)/pList[vAB-1])/((double) steps);
		}
		}

/*CIs[listNum*3-3] = pList[1];
CIs[listNum*3-2] = pList[99];
CIs[listNum*3-1] = listSum;
*/	}


    /**
	TODO: enable CI calculation for indels
     */
    private static int[][] getIndelContingencyTable(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, final Map<String, PerReadAlleleLikelihoodMap> stratifiedLikelihoodMap) {
        final double INDEL_LIKELIHOOD_THRESH = 0.3;

        if (stratifiedLikelihoodMap == null)
            return null;
        
        int[][] table = new int[2][2];

        for ( String sample : stratifiedContexts.keySet() ) {
            final PerReadAlleleLikelihoodMap indelLikelihoodMap = stratifiedLikelihoodMap.get(sample);
            final AlignmentContext context = stratifiedContexts.get(sample);
            if ( context == null || indelLikelihoodMap == null ) 
                continue;

            ReadBackedPileup pileup = context.getBasePileup();
            if (pileup == null)
                 continue;

            for (final PileupElement p: pileup) {
                if ( p.getRead().getMappingQuality() < 20)
                    continue;
                if (indelLikelihoodMap.containsPileupElement(p)) {
                    // to classify a pileup element as ref or alt, we look at the likelihood associated with the allele associated to this element.
                    // A pileup element then has a list of pairs of form (Allele, likelihood of this allele).
                    // To classify a pileup element as Ref or Alt, we look at the likelihood of corresponding alleles.
                    // If likelihood of ref allele > highest likelihood of all alt alleles  + epsilon, then this pileup element is "ref"
                    // otherwise  if highest alt allele likelihood is > ref likelihood + epsilon, then this pileup element it "alt"
                    // retrieve likelihood information corresponding to this read
                    Map<Allele,Double> el = indelLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(p);
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
