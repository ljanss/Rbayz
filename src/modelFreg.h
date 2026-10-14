//
//  BayzR -- modelFreg.h
//  Computational class defining methods to work on one covariate.
//  - the parameter vector has length 1, and no levelNames.
//  - no hyper parameters - sampleHpar is empty.
//  - also defines sample(). So far only one model-object that uses one
//    covariate, so it is combining all de/correction, lhs/rhs statistics
//    and sample code.
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelFreg_h
#define modelFreg_h

#include <Rcpp.h>
#include <cmath>
#include "modelCoeff.h"
#include "dataCovar.h"
#include "parsedModelTerm.h"

class modelFreg : public modelCoeff {
   
public:
   
   modelFreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod),
           C(modeldescr.variableObjects[0]),
           comp(C)
   {
      if(modeldescr.variableTypes[0] != 3)
         throw generalRbayzError("Not a numeric vector input: "+modeldescr.variableNames[0]);
      par = new parVector(modeldescr, 0.0l);
      par->traced=1;
   }
   
   ~modelFreg() {
      delete par;
   }
   
   void sample() {
      resid.add(fit);
      comp.collectLhsRhs(resid, residPrec, lhs, rhs);
      par->val[0] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
      comp.fillFit(fit, par->val[0]);
      resid.sub(fit);
   }

   void sampleHpars() {}

   void restart() {}

   dataCovar C;
   dataAlgoCovar comp;
   double lhs, rhs;
   
};

#endif /* modelFreg_h */
