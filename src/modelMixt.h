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
      varmodel = new scaleVarStr(modeldescr, this->par);
      p_incl = 0.05;
      p_excl = 0.20;
   }

   ~modelMixt() {
      delete varmodel;
      delete pi;
   }
   
   void sample() {
      // Start with relatively simple Include-Move-Exclude strategy.
      double inf = std::numeric_limits<double>::infinity();
      double lhs, rhs;
      for(size_t k=0; k < M.ncol; k++) {
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
   indepVarStr* varmodel;
   double p_incl, p_excl;

};

#endif /* modelMixt_h */
