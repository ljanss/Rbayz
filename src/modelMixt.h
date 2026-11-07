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
      mxi = new parVector(modeldescr, 0.0l, *par, "mxi");
      Rbayz::parList.push_back(&mxi);
      // The mixture model can have difficulty starting when initial variance is too
      // large (and the default start value of 1 is often far too large). This causes no
      // acceptances / switches of beta's to the '1' class, all betas stay 0, and fit stays 0.
      // For now I'm trying to get around this problem by computing a sensible initial starting
      // variance below, and on top of that the sampleScale() method of the variance class
      // will make up to 20 shrinkage adjustments (and then gives up with an error).
      // If all this is not enough, the alternative would be an initial sampling of beta's
      // and mixture indicators to get a non-zero fit from the start.
      double initvar = 0.5 * M.sumColVar / double(M.ncol);
      varmodel = new scaleVarStr(modeldescr, this->par, initvar);
      // parameters for MH sampler (under development)
      p_incl = 0.05;
      p_excl = 0.20;
   }

   ~modelMixt() {
      delete varmodel;
      delete pi;
   }
   
   void sample() {
      // Start with Daniel's full gibbs implementation for a 2-mixture model.
      // The set-up for the mixture model is now that one class handles two sets of parameters:
      // the mixture indicators (mxi) and the 'betas' (par vector). And similarly in the hyper
      // paramters it handles variance and pi's.
      // Formally updating mxi and beta is two different tasks doing mxi|beta, and then
      // beta|mxi, but done here per-column in a single loop. The first part is the mxi|beta
      // update. Here if current mxi is zero, the generated new_beta is formally the current
      // beta (it could have been generated the previous round and stored, but it's
      // enough to make it just before it is needed). The second part in the loop then
      // updates beta|mxi, which considers all mxi that are 1 (either just set to 1 or
      // already 1 before). This second loop can skip all cases where mxi is 0, as it implies
      // beta is zero.
      // Note: this now uses a model with scale, so it is fitting scale*sum(b_i x_i), and
      // the beta samples are on N(0,1) scale. So there is some juggling with this scale:
      // the betaUpdate and resid_correct/decorrect functions don't know about this scale,
      // and it must be explicitly added; the beta stored, though, is on the N(0,1) scale;
      // the collectLhsRhs function also don't include the scale when making their statistics,
      // and scale needs to be added afterwards as an extra scale^2 in lhs, and an extra scale
      // in rhs (considering that the covariate used is scale*x_i). This simplifies in
      // rhs/lhs to an extra divition by scale, and in sqrt(1.0/lhs) can add an extra division
      // by scale outside the sqrt(). For fit-functions there are versions that do take the scale
      // into account and can work with the N(0,1) beta's. Note: fillFit with scale skips
      // zero beta's but not other fit-functions!
//      double inf = std::numeric_limits<double>::infinity();
      double lhs, rhs;
      double sse0, sse1;
      double beta_new, temp;
      double log_prob_ratio;
      double log_pi_ratio = log(pi->val[0]/pi->val[1]);
      double scale = sqrt(varmodel->par->val[0]);
      for(size_t k=0; k < M.ncol; k++) {
         // Part 1: updating mxi[i]|beta[i]. 
         // Note: sse_eval_betaUpdate() expects beta_diff = old - new value, so
         // if old=0, new=beta_new, pass -beta_new; if old is par-val[k] and
         // new=0, pass par->val[k].
         if (mxi->val[k] == 0.0l) {
            sse0 = collect_sse();
            beta_new = R::rnorm(0, 1);
            sse1 = comp.sse_eval_betaUpdate(resid, residPrec, fit, -beta_new, k);
         }
         else {
            sse0 = comp.sse_eval_betaUpdate(resid, residPrec, fit, par->val[k], k);
            sse1 = collect_sse();
         }
         log_prob_ratio = -0.5 * (sse1 - sse0) + log_pi_ratio;
         temp = R::runif(0, 1);
         if (log(temp/(1-temp)) < log_prob_ratio) { // set mxi[k] to 0.0 (if it was 1)
            if(mxi->val[k] == 1.0l) {
               mxi->val[k] = 0.0l;
               comp.resid_decorrect(resid, par->val[k], k);
               par->val[k] = 0.0l;
            }
         }
         else { // set mxi[k] to 1.0 )(if it was zero)
            if(mxi->val[k] == 0.0l) {
               mxi->val[k] = 1.0l;
               par->val[k] = beta_new;
               comp.resid_correct(resid, par->val[k], k);
            }
         }
         // Part 2: update beta[i]|mxi[i] (only if mxi[i] is 1)
         if (mxi->val[k] == 1.0l) {
            comp.resid_decorrect(resid, (temp=scale*par->val[k]), k);
            comp.collectLhsRhs(resid, residPrec, lhs, rhs, k);
            lhs += 1.0l;
            par->val[k] = R::rnorm( (rhs/(lhs*scale)), (sqrt(1.0/lhs)/scale));
            comp.resid_correct(resid, (temp=scale*par->val[k]), k);
         }
      }
      comp.fillFit(fit, par->Values, scale);
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
