//
// BayzR -- dataAlgoIntercept.h
//
// Algorithms to collect statistics for updating intercept (mean).
//
// Created by Luc Janss on 13/10/2026.

#ifndef dataAlgoIntercept_h
#define dataAlgoIntercept_h

#include "simpleVector.h"
#include "dataAlgorithm.h"

class dataAlgoIntercept : public dataAlgorithm {
   
public:

   dataAlgoIntercept() = default;

   ~dataAlgoIntercept() = default;

   void fillFit(simpleDblVec& fit, simpleDblVector& beta) {
      for (size_t obs=0; obs < fit.nelem; obs++)
        fit[obs] = beta->val[0];
   }

   void fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale) {
      for (size_t obs=0; obs < fit.nelem; obs++)
        fit[obs] = beta->val[0] * scale;
   }

   void resid_correct(simpleDblVec& resid, simpleDblVector& beta) {
      for (size_t obs=0; obs < resid.nelem; obs++)
        resid[obs] -= beta->val[0];
   }

   void resid_decorrect(simpleDblVec& resid, simpleDblVector& beta) {
      for (size_t obs=0; obs < resid.nelem; obs++)
        resid[obs] += beta->val[0];
   }

   void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs) {
      lhs[0] = 0.0;
      rhs[0] = 0.0;
      for (size_t obs=0; obs < resid.nelem; obs++) {
         lhs[0] += residprec[obs];
         rhs[0] += resid[obs]*residprec[obs];
      }
   }

};

#endif /* dataAlgoIntercept_h */