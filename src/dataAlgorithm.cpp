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

void dataAlgoIntercept::resid_correct(simpleDblVec& resid, simpleDblVector& beta) {
   for (size_t obs=0; obs < resid.nelem; obs++)
     resid[obs] -= beta->val[0];
}

void dataAlgoIntercept::resid_decorrect(simpleDblVec& resid, simpleDblVector& beta) {
   for (size_t obs=0; obs < resid.nelem; obs++)
     resid[obs] += beta->val[0];
}

void dataAlgoIntercept::collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs) {
   lhs[0] = 0.0;
   rhs[0] = 0.0;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      lhs[0] += residprec[obs];
      rhs[0] += resid[obs]*residprec[obs];
   }
}

