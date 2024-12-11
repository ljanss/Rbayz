//
//  Rbayz --- modelMixt.h
//  Object to modify mixture class indicator for MIXT variance structure.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMixt_h
#define modelMixt_h

#include <Rcpp.h>
#include "modelBase.h"

class modelMixt : public modelBase {

public:

   modelMixt(parsedModelTerm & modeldescr, modelRreg * rrmod)
         : modelBase()
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
