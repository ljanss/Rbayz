// BayzR -- modelRanfi.h
// Model class for RANdom Factor with Indep variance structures.
// Derives from modelFactor and has pointer to IndepVarStr object to model the variance.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanfi_h
#define modelRanfi_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "indepVarStr.h"
//#include <unistd.h>

class modelRanfi : public modelFactor {

public:

   modelRanfi(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod, true) {   // the true is for collapseInteractions
   }

   ~modelRanfi() {                 // the varmodels get allocated in the derived classes, but
      if (varmodel != 0)           // by putting destructor here it does not have to repeat
         delete varmodel;          // the same destructor in each derived class
   }

   void sample() {
      // sample() method for random effects with indep var-structure: the data
      // corrections and LHS and RHS can be made using the methods from modelFactor.
      // To update parameters as random effects varmodel->weights[k] are added in LHS.
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=0; k<par->nelem; k++) {
         lhs[k] += varmodel->weights[k];
         par->val[k] = R::rnorm((rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
      }
      resid_correct();
   }

   void sampleHpars() {
      varmodel->sample();
   }

   void restart() {
      varmodel->restart();
   }

   indepVarStr* varmodel=0;

};

// Implementation for IDEN variance structure
// Note: all implementations with different INDEP variance structures can use methods from
// the modelRanfi parent class, there is only need for allocating the right varmodel object.
// Futher extensions with indep structures can include diag, weighted, loglin, mixed?

class modelRanfi_iden : public modelRanfi {
public:
   modelRanfi_iden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRanfi(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};


#endif /* modelRanfi_h */
