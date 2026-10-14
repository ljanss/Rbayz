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

// Implementations of dataAlgoFactor methods
// -----------------------------------------

void dataAlgoFactor::fillFit(simpleDblVector& fit, simpleDblVector& beta) {
   for (size_t obs=0; obs < Fdata.nelem; obs++)
     fit[obs] = beta->val[Fdata.data[obs]];
}

void dataAlgoFactor::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < Fdata.nelem; obs++)
     fit[obs] = beta->val[Fdata.data[obs]] * scale;
}

void dataAlgoFactor::collectLhsRhs(simpleDblVector& resid, simpleDblVector& residprec,
                                   simpleDblVector& lhs, simpleDblVector& rhs) {
   rhs.fill(0.0);
   lhs.fill(0.0);
   for (size_t obs=0; obs < Fdata.nelem; obs++) {
      size_t k = Fdata.data[obs];
      rhs.data[k] += residprec[obs] * resid[obs];
      lhs.data[k] += residprec[obs];
   }
}

// Implementations of dataAlgoCovar methods
// -----------------------------------------

void dataAlgoCovar::fillFit(simpleDblVector& fit, double& beta) {
   for (size_t obs=0; obs < Cdata.nelem; obs++)
     fit[obs] = beta * Cdata.data[obs];
}

void dataAlgoCovar::fillFit(simpleDblVector& fit, double& beta, double scale) {
   for (size_t obs=0; obs < Cdata.nelem; obs++)
     fit[obs] = beta * Cdata.data[obs] * scale;
}

void dataAlgoCovar::collectLhsRhs(simpleDblVector& resid, simpleDblVector& residprec, double& lhs, double& rhs) {
   lhs = 0.0;
   rhs = 0.0;
   for (size_t obs=0; obs < Cdata.nelem; obs++) {
      lhs += residprec[obs] * Cdata.data[obs] * Cdata.data[obs];
      rhs += residprec[obs] * resid[obs] * Cdata.data[obs];
   }
}

// Implementations of nested regressions ... to be done


