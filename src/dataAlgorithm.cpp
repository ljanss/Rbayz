//
// BayzR -- dataAlgorithm.cpp
//
// Implementations of dataAlgorithm and derived classes methods.
//
// Created by Luc Janss on 13/10/2026.

#include "dataAlgoIntercept.h"

// Implementations of dataAlgoIntercept methods
// --------------------------------------------

void dataAlgoIntercept::fillFit(simpleDblVec& fit, simpleDblVector& beta) {
   for (size_t obs=0; obs < fit.nelem; obs++)
     fit[obs] = beta->val[0];
}

void dataAlgoIntercept::fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < fit.nelem; obs++)
     fit[obs] = beta->val[0] * scale;
}

void dataAlgoIntercept::collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, double& lhs_sc, double& rhs_sc) {
   lhs_sc = 0.0;
   rhs_sc = 0.0;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      lhs_sc += residprec[obs];
      rhs_sc += resid[obs]*residprec[obs];
   }
}

