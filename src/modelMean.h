//
//  BayzR --- modelMean.h
//
// Computational class for fitting mean (intercept)
// - this one has no data!
// - par is size 1 to hold mean
// - hpar is not used
// Note: the mean-term goes through the parsing because there is a "1" inserted in the
// model-terms, but it means the name set in the base-class constructor is "1" and needs
// to be changed.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMean_h
#define modelMean_h

#include <Rcpp.h>
#include "modelCoeff.h"

class modelMean : public modelCoeff {
   
public:

   modelMean(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod), comp() {
      par = new parVector(modeldescr, 0.0l);
      par->Name="mean";       // a little repair because the
      par->Labels[0]="mean";  // constructor has not set it right
      par->traced=1;
   }

   ~modelMean() {
      delete par;
   }

   void sample() {
      resid.add(fit);
      comp.collectLhsRhs(resid, residPrec, lhs, rhs);
      par->val[0] = R::rnorm((rhs[0]/lhs[0]), sqrt(1.0/lhs[0]));
      comp.fillFit(fit, par);
      resid.sub(fit);
   }

   void sampleHpars() {}

   void restart() {}

   dataAlgoIntercept comp;

private:
};

#endif /* modelMean_h */
