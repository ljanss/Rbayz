//
//  Rbayz --- modelMixt.h
//  Object to modify mixture class indicator for MIXT variance structure.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMixt_h
#define modelMixt_h

#include <Rcpp.h>
#include "modelMatrix.h"

class modelMixt : public modelBase {

public:

   modelMixt(parsedModelTerm & modeldescr, modelRreg * rrmod)
         : modelbase()
   {
      par = regcoeff;

   }

   ~modelMixt() {
   }
   
   void sample() {
      update_regressions();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]/M->weights[k];
      hpar[0] = gprior.samplevar(ssq, M->ncol);
   }

   void sampleHpars() {

   }

   void restart() {}

};

#endif /* modelMixt_h */
