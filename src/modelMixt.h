//
// Rbayz --- modelMixt.h
// Modelling of Variable Selection / mixture prior models, now implemented with its
// own model-term 'mx'.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMixt_h
#define modelMixt_h

#include <Rcpp.h>
#include "modelBase.h"

class modelMixt : public modelMatrix {

public:

   modelMixt(parsedModelTerm & modeldescr, modelResp* rmod)
         : modelMatrix(modeldescr, rmod), pi()
   {
      // These problems should be caught already in model parsing ....
      if (! ( modeldescr.allOptions["vars"].isgiven && modeldescr.allOptions["counts"].isgiven) ) {
         throw generalRbayzError("Options 'vars' and 'counts' are required for mx() model term");
      }
      vars = modeldescr.allOptions["vars"].valnumb;
      counts = modeldescr.allOptions["counts"].valnumb;
      if (vars.size() != counts.size()) {
         throw generalRbayzError("The number of 'vars' must match the number of 'counts' for mx() model term");
      }
      Nmixtures = vars.size();
      // For the moment, only work with 2 mixtures where var[0]=0 and var[1]=1
      if(Nmixtures != 2) {
         throw generalRbayzError("Rbayz is not yet ready to work with more than 2 mixtures ...");
      }
      vars[0] = 0.0l;
      vars[1] = 1.0l;
      if (counts[0] <1 || counts[1] < 1) {
         throw generalRbayzError("Each 'count' for mx() model term must be at least 1");
      }
      std::vector<std::string> piNames(2);
      piNames[0] = "pi0." + modeldescr.variableString;
      piNames[1] = "pi1." + modeldescr.variableString;
      pi = new parVector(modeldescr, 0.5l, piNames);   // extra parameter that is here manually
      Rbayz::parList.push_back(&pi);                   // added in the global parList!
      double sum_counts = counts[0] + counts[1];
      pi->val[0] = counts[0] / sum_counts;
      pi->val[1] = counts[1] / sum_counts;
      // little dount here if mxi vector could best be integer or double;
      // when averaged it gets double, but as integer (or simpleSizeVector) or
      // can be used directly as indexes to indicate mixture group.
      // Or it could be both - then use prepForOutput() method to copy from
      // internal Size vector or 'public' double par vector.
      mxi = new parVector(modeldescr, 0.0, par, "mxi");
      Rbayz::parList.push_back(&mxi);
      // Compute a sensible initial variance for the scaleVarStr class, note it is variance!
      double initvar = 0.5 * M.sumColVar / double(M.ncol);
      varmodel = new scaleVarStr(modeldescr, this->par, initvar);
      p_incl = 0.05;
      p_excl = 0.20;
      // This mixture model may have difficulty starting when all initial
      // mxi are zero, and initial beta variance is too big (e.g. default start 1).
      // Option 1: initialize some random mxi at 1, even though effect samples
      // are difficult to make realistic? This should trigger variance updates to
      // make variance much smaller.
      // Option 2: try to set some sensible initial variance. I could maybe manage
      // that in the scale-sampling algorithm.

   }

   ~modelMixt() {
      delete varmodel;
      delete pi;
   }
   
   void sample() {
      // Start with Daniel's full gibbs implementation, this is
      // two mixture model with spike at zero.
      double inf = std::numeric_limits<double>::infinity();
      double lhs, rhs;
      for(size_t k=0; k < M.ncol; k++) {
         if (mxi[k] == 0.0l) {

         }
         else {

         }
         comp.resid_decorrect(resid, par->val[k], k);
         if(varmodel->weights[k]==inf) {
            par->val[k]=0;
         }
         else {
            comp.collectLhsRhs(resid, residPrec, lhs, rhs, k);
            lhs += varmodel->weights[k];
            par->val[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         }
         comp.resid_correct(resid, par->val[k], k);
      }
      comp.fillFit(fit, par->Values);
   }

   void sampleHpars() {
      // note getFitScaleStats and resid_fit_scaleUpdate are parts of modelCoeff
      // class, and not part of the dataAlgo 'comp' object.
      double oldscale = sqrt(varmodel->par->val[0]);
      double lhs, rhs;
      getFitScaleStats(lhs, rhs);
      double newscale = varmodel->sampleScale(lhs, rhs);
      resid_fit_scaleUpdate(oldscale,newscale);
      // Also update pi's. There is now not a dedicated class that models pi's,
      // so it's just added here in the modelMixt object.

   }

   void restart() {}

   size_t Nmixtures;
   std::vector<double> vars;
   std::vector<double> counts;
   parVector* pi;
   parVector* mxi;
   indepVarStr* varmodel;
   double p_incl, p_excl;

};

#endif /* modelMixt_h */
