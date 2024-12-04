//
//  BayzR --- modelRreg.hpp
//
//  Computational class to model random regressions on a matrix of covariates.
//  This uses methods from modelMatrix, only some parameter names need to be set.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRreg_h
#define modelRreg_h

#include <Rcpp.h>
#include "modelMatrix.h"
#include "indepVarStr.h"
#include "dataMatrix.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelMatrix(modeldescr, rmod)   {
      if(checkOptions(modeldescr.options, "V prior pvals")>0) {
         throw(generalRbayzError("ERROR: unrecognized option(s) in "+modeldescr.shortModelTerm));
      }
   }

   ~modelRreg() {
      delete varmodel;
   }
   
   // sample() now updated to handle zero variance (inf weight) from the variance model:
   // - decorrect() runs to correct for old estimate, but skips if old one was already zero
   // - for inf weight, regcoeff is set to zero, else it will go through the sampling procedure
   // - correct() runs in case regcoeff != zero, but skips if it is zero
   void sample() {
      double inf = std::numeric_limits<double>::infinity();
      for(size_t k=0; k < M->ncol; k++) {
         resid_decorrect(k);
         if(varmodel->weights[k]==inf) {
            par->val[k]=0;
         }
         else {
            collect_lhs_rhs(k);   // update lhs and rhs variables
            lhs += varmodel->weights[k];
            par->val[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         }
         resid_correct(k);
      }
      // need some thinking how to store sample info in file; likely every "saved" cycle.
      // main runs prepForOutput, and on a parVector main runs collectStats at the save intervals,
      // or it needs a new mechanism to switch on saving samples from output (which can be generic feature).
   }

   void sampleHpars() {
      varmodel->sample();
   }

   indepVarStr* varmodel;

};

// Here can start working on defining different variance structures for modelRreg
// For the moment, Rreg is only designed to accept variance models in the "indepVarStr" class, but this
// can include LASSO, BVS, DIAG / weighted, log-linear ...
class modelRregIden : public modelRreg {
public:
   modelRregIden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};

class modelRregDiag : public modelRreg {
public:
   modelRregDiag(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new diagVarStr(pmdescr, this->par);
   }
};

class modelRregMixt : public modelRreg {
public:
   modelRregMixt(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new mixtVarStr(pmdescr, this->par);
   }
};

#endif /* modelRreg */
