//
//  Rbayz --- modelMixt.h
//  Modelling of Variable Selection / mixture prior models.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMixt_h
#define modelMixt_h

#include <Rcpp.h>
#include "modelBase.h"

class modelMixt : public modelMatrix {

public:

   modelMixt(parsedModelTerm & modeldescr, modelRreg * rrmod)
         : modelMatrix(modeldescr, rrmod)
   {

   }

   ~modelMixt() {
   }
   
   void sample() {
      // update_regressions();
      // update hyper-par (variance) using SSQ of random effects
   }

   void sampleHpars() {

   }

   void restart() {}

};

#endif /* modelMixt_h */
