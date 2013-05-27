/* Base class extensing variant walkers, to be derived from Clojure.
 */

package bcbio.variation;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
//import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
//import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
//import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
//import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;
import java.io.PrintStream;

public abstract class BaseVariantWalker extends RodWalker {
  @Output
  public String out;

  @ArgumentCollection
  public StandardVariantContextInputArgumentCollection invrns = new StandardVariantContextInputArgumentCollection();

//   @Input(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", required=true)
//   public List<RodBinding<VariantContext>> evals;

//   @Override
//   public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
// 	return 0;
//   }

//   @Override
//   public Integer reduceInit () {
//     return 0;
//   }

//   @Override
//   public Integer reduce(Integer value, Integer sum) {
//     return 0;
//   }

//   @Override
//   public void onTraversalDone (Integer result) {
//   }
}





