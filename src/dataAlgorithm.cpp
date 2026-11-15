//
// BayzR -- dataAlgorithm.cpp
//
// Implementations of dataAlgorithm and derived classes methods.
//
// Created by Luc Janss on 13/10/2026.

#include "dataAlgorithm.h"

// Utility functions
// -----------------
// buildObsIndex used in dataAlgoMatrix and a helper function for searching.

// Comparison function to search for a string in a vector of pair<string, number>
inline bool compString2Pair(const std::pair<std::string, size_t> & p, const std::string & s)
{
    return p.first < s;
}

void builObsIndex(simpleSizeVector& obsIndex, dataFactor& F, labeledMatrix& M) {
   int errors=0;
   // Build a sorted list of the matrix rownames paired with matrix entry-rows
   std::vector< std::pair<std::string, size_t> > matLabelsSorted;
   for(size_t i=0; i< M.rownames.size(); i++)
      matLabelsSorted.push_back(std::make_pair(M.rownames[i], i));
   std::sort(matLabelsSorted.begin(),matLabelsSorted.end()); // shoud be enough
   // Make the index that links observations to matrix rows
   std::string s;
   std::vector< std::pair<std::string, size_t> >::iterator it;
   for(size_t i=0; i<F.nelem; i++) {
      s = F.labels[F.data[i]];
      if ( (it = std::lower_bound(matLabelsSorted.begin(), matLabelsSorted.end(), s, compString2Pair)) != matLabelsSorted.end()) {
         obsIndex[i]=it->second;
      }
      else {
         // error [ToDo] can be nicer in messages list
         Rcpp::Rcout << "ID " << s << " is not available in similarity/kernel matrix\n";
         errors++;
      }
   }
   if(errors)
      throw (generalRbayzError("Some IDs not found in similarity/kernel matrix"));
}

// Implementations of dataAlgoIntercept methods
// --------------------------------------------

void dataAlgoIntercept::fillFit(simpleDblVector& fit, simpleDblVector& beta) {
   for (size_t obs=0; obs < fit.nelem; obs++)
     fit[obs] = beta.data[0];
}

void dataAlgoIntercept::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < fit.nelem; obs++)
     fit[obs] = beta.data[0] * scale;
}

void dataAlgoIntercept::collectLhsRhs(simpleDblVector& resid, simpleDblVector& residprec, double& lhs_sc, double& rhs_sc) {
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

void dataAlgoCovar::fillFit(simpleDblVector& fit, simpleDblVector& beta) {
   for (size_t obs=0; obs < Cdata.nelem; obs++)
     fit[obs] = beta.data[0] * Cdata.data[obs];
}

void dataAlgoCovar::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   for (size_t obs=0; obs < Cdata.nelem; obs++)
     fit[obs] = beta.data[0] * Cdata.data[obs] * scale;
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

// pre-computing the x'x diagonals for simple cases. residdiag is an optional diagonal
// matrix from the residual variance, it can be nullptr and the interface takes pointer
// instead of reference to allow for this.
// ***************
// **** OBS: not yet correct!
// **** This should also use obsindx to only use Mdata rows that are matched to residuals.
// ***************
void dataAlgoMatrix::comp_xx_diag(simpleDblVector& xxdiag, const simpleDblVector* residdiag) {
   double * colptr;
   if(residdiag == nullptr) {
      for(size_t k=0; k < Mdata.ncol; k++) {
         colptr = Mdata.data[k];
         xxdiag.data[k] = 0.0;
         for(size_t obs=0; obs < Mdata.nrow; obs++)
            xxdiag.data[k] += colptr[obs] * colptr[obs];
      }
   }
   else {
      for(size_t k=0; k < Mdata.ncol; k++) {
         colptr = Mdata.data[k];
         xxdiag.data[k] = 0.0;
         for(size_t obs=0; obs < Mdata.nrow; obs++)
            xxdiag.data[k] += colptr[obs] * colptr[obs] * residdiag->data[obs];
      }
   }
}


void dataAlgoMatrix::fillFit(simpleDblVector& fit, simpleDblVector& beta) {
   double * colptr;
   for (size_t obs=0; obs < fit.nelem; obs++) fit.data[obs] = 0.0l;
   for(size_t k=0; k < Mdata.ncol; k++) {
      if (beta.data[k] == 0.0l) continue;
      colptr = Mdata.data[k];
      for (size_t obs=0; obs < fit.nelem; obs++)
         fit.data[obs] += beta.data[k] * colptr[obsIndex[obs]];
   }
}

// This fillFit() with scale skips zero beta's for efficiency. The case with (many) zero
// beta's occurs on mixture models which are so far the only class using this fillFit().
void dataAlgoMatrix::fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale) {
   double * colptr;
   for (size_t obs=0; obs < fit.nelem; obs++) fit.data[obs] = 0.0l;
   for(size_t k=0; k < Mdata.ncol; k++) {
      if (beta.data[k] == 0.0l) continue;
      colptr = Mdata.data[k];
      for (size_t obs=0; obs < fit.nelem; obs++)
         fit.data[obs] += beta.data[k] * colptr[obsIndex[obs]] * scale;
   }
}

void dataAlgoMatrix::accumFit(simpleDblVector& fit, double& beta, size_t& col) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < fit.nelem; obs++)
      fit.data[obs] += beta * colptr[obsIndex[obs]];
}

void dataAlgoMatrix::accumFit(simpleDblVector& fit, double& beta, size_t& col, double scale) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < fit.nelem; obs++)
      fit.data[obs] += beta * colptr[obsIndex[obs]] * scale;
}

void dataAlgoMatrix::resid_correct(simpleDblVector& resid, double beta, size_t col) {
   if(beta==0.0l) return;
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++) {
      resid.data[obs] -= beta * colptr[obsIndex[obs]];
   }
}

void dataAlgoMatrix::resid_decorrect(simpleDblVector& resid, double beta, size_t col) {
   if(beta==0.0l) return;
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++)
      resid.data[obs] += beta * colptr[obsIndex[obs]];
}

// adjusting residuals and fit for a proposed change in beta (beta_diff)
void dataAlgoMatrix::resid_fit_betaUpdate(simpleDblVector& resid, simpleDblVector& fit, double beta_diff, size_t col) {
   double * colptr = Mdata.data[col];
   for (size_t obs=0; obs < resid.nelem; obs++) {
      resid.data[obs] += beta_diff * colptr[obsIndex[obs]];
      fit.data[obs] -= beta_diff * colptr[obsIndex[obs]];
   }
}

// A version that has merged the residual and fit adjustments AND collects LHS and RHS
// to update beta. The beta coming in is the old one that is both removed (+) from residual
// and from fit (-). After updating beta, this is combined with resid_fit_betaUpdate to
// simultaneously update resid and fit again with the new beta.
void dataAlgoMatrix::merged_adjust_getLhsRhs(simpleDblVector& resid, simpleDblVector& residPrec,
         simpleDblVector& fit, double beta, double& lhs, double& rhs, size_t col) {
   lhs = 0.0l; rhs = 0.0l;
   double * colptr = Mdata.data[col];
   double x, temp1;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      x = colptr[obsIndex[obs]];
      resid.data[obs] += beta * x;
      fit.data[obs] -= beta * x;
      temp1 = x * residPrec[obs];
      rhs += temp1 * resid.data[obs];
      lhs += temp1 * x;
   }
}

// merge_adjust_getLhsRhs versions specialized for IDEN and DIAG structures.
// Both use the pre-computed x'x (x'Wx for DIAG) LHS diagonals, only need to multiple
// with the inverse scalar variance. For the RHS computation, IDEN skips diagonal weights,
// but DIAG still needs (now using the prepared residPrec as in other computations).
// OBS: these methods depend on the calling function (sample()) storing the current
// residual variance in the member residScalarVar.
void dataAlgoMatrix::merged_adjust_getLhsRhs_iden(simpleDblVector& resid, simpleDblVector& residPrec,
         simpleDblVector& fit, double beta, double& lhs, double& rhs, size_t col) {
   lhs = xx_diag[col] / residScalarVar;
   rhs = 0.0l;
   double * colptr = Mdata.data[col];
   double x;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      x = colptr[obsIndex[obs]];
      resid.data[obs] += beta * x;
      fit.data[obs] -= beta * x;
      rhs += x * resid.data[obs];
   }
   rhs /= residScalarVar;
}

void dataAlgoMatrix::merged_adjust_getLhsRhs_diag(simpleDblVector& resid, simpleDblVector& residPrec,
         simpleDblVector& fit, double beta, double& lhs, double& rhs, size_t col) {
   lhs = xx_diag[col] / residScalarVar;
   rhs = 0.0l;
   double * colptr = Mdata.data[col];
   double x;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      x = colptr[obsIndex[obs]];
      resid.data[obs] += beta * x;
      fit.data[obs] -= beta * x;
      rhs += x * residPrec[obs] * resid.data[obs];
   }
   // note rhs uses here residPrec that already includes the inverse variance
}

void dataAlgoMatrix::sse_collectCurrentAndBetaUpdate(simpleDblVector& resid, simpleDblVector& residPrec,
         double beta_diff, double& sse_current, double& sse_betaUpdate, size_t col) {
   double * colptr = Mdata.data[col];
   double temp0, temp1;
   sse_current = 0.0l;
   sse_betaUpdate = 0.0l;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      temp0 = resid.data[obs];
      temp1 = temp0 + beta_diff * colptr[obsIndex[obs]];
      sse_current += temp0 * temp0 * residPrec[obs];
      sse_betaUpdate += temp1 * temp1 * residPrec[obs];
   }
}

// Evaluate SSE for a (proposed) change in beta. It uses the same sign convention as
// resid_fit_betaUpdate where the beta_diff to pass is (old-new) value.
double dataAlgoMatrix::sse_eval_betaUpdate(simpleDblVector& resid, simpleDblVector& residPrec,
            simpleDblVector& fit,  double beta_diff, size_t col)
{
   double * colptr = Mdata.data[col];
   double sse = 0.0l, temp=0.0l;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      temp = resid.data[obs] + beta_diff * colptr[obsIndex[obs]];
      sse += temp * temp * residPrec[obs];
   }
   return sse;
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

void dataAlgoMatrix::collectLhsRhs_withAccumFit(simpleDblVector& resid, simpleDblVector& residPrec,
         simpleDblVector& fit, double& lhs, double& rhs, size_t col) {
   size_t matrixrow;
   lhs = 0.0l; rhs = 0.0l;
   double * colptr = Mdata.data[col];
   double temp1;
   for (size_t obs=0; obs < resid.nelem; obs++) {
      matrixrow = obsIndex[obs];
      temp1 = colptr[matrixrow] * residPrec[obs];
      rhs += temp1 * (resid.data[obs] - fit.data[obs]);
      lhs += temp1 * colptr[matrixrow];
   }
}
