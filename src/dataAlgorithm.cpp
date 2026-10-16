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
     fit[obs] = beta.data[0];
}

void dataAlgoIntercept::fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < fit.nelem; obs++)
     fit[obs] = beta.data[0] * scale;
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
     fit[obs] = beta.data[Fdata.data[obs]];
}

void dataAlgoFactor::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < Fdata.nelem; obs++)
     fit[obs] = beta.data[Fdata.data[obs]] * scale;
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

// Implementations of dataAlgoMatrix methods
// -----------------------------------------

void dataAlgoMatrix::fillFit(simpleDblVector& fit, simpleDblVector& beta) {
   double * colptr;
   for (size_t obs=0; obs < fit.nelem; obs++) fit.data[obs] = 0.0l;
   for(size_t k=0; k < Mdata.ncol; k++) {
      colptr = Mdata.data[k];
      for (size_t obs=0; obs < fit.nelem; obs++)
         fit.data[obs] += beta.data[k] * colptr[Fdata.data[obs]];
   }
}

void dataAlgoMatrix::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   double * colptr;
   for (size_t obs=0; obs < fit.nelem; obs++) fit.data[obs] = 0.0l;
   for(size_t k=0; k < Mdata.ncol; k++) {
      colptr = Mdata.data[k];
      for (size_t obs=0; obs < fit.nelem; obs++)
         fit.data[obs] += beta.data[k] * colptr[Fdata.data[obs]] * scale;
   }
}

void dataAlgoMatrix::accumFit(simpleDblVector& fit, double& beta, size_t& col) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < Fdata.nelem; obs++)
      fit.data[obs] += beta * colptr[Fdata.data[obs]];
}

void dataAlgoMatrix::accumFit(simpleDblVector& fit, double& beta, size_t& col, double scale) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < fit.nelem; obs++)
      fit.data[obs] += beta * colptr[Fdata.data[obs]] * scale;
}

void dataAlgoMatrix::resid_correct(simpleDblVector& resid, double& beta, size_t col) {
   if(beta==0.0l) return;
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++) {
      resid.data[obs] -= beta * colptr[Fdata.data[obs]];
   }
}
void dataAlgoMatrix::resid_decorrect(simpleDblVector& resid, double& beta, size_t col) {
   if(beta==0.0l) return;
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++)
      resid.data[obs] += beta * colptr[Fdata.data[obs]];
}
void dataAlgoMatrix::resid_fit_betaUpdate(simpleDblVector& resid, double beta_diff, size_t col) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++) {
      resid.data[obs] += beta_diff * colptr[obsIndex[obs]];
      fit.data[obs] -= beta_diff * colptr[obsIndex[obs]];
   }
}

void dataAlgoMatrix::collectLhsRhs(simpleDblVector& resid, simpleDblVector& residPrec,
         double& lhs, double& rhs, size_t col) {
   size_t matrixrow;
   lhs = 0.0l; rhs = 0.0l;
   double * colptr = Mdata.data[col];
   double temp1;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      matrixrow = obsIndex[obs];
      temp1 = colptr[matrixrow] * residPrec[obs];
      rhs += temp1 * resid.data[obs];
      lhs += temp1 * colptr[matrixrow];
   }
}

