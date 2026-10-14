//  R/bayz
//  modelMatrix.h
//
//  Defines the computational methods when design matrix is a (real) matrix of covariates.
//  Also needs a factor to make the link to data, but is not derived from modelFactor because
//  1) modelFactor would allocate a par-vector at the wrong size; 2) modelFactor can also handle
//  interactions, but that is not supported in modelMatrix (it would complicate linking to the data
//  because the interaction coded in modelFactor does not exist as a single variable in the data).
//
//  Created by Luc Janss on 30/08/2019.
//

#ifndef modelMatrix_h
#define modelMatrix_h

#include <Rcpp.h>
#include <cmath>
#include "dataFactor.h"
#include "dataMatrix.h"
#include "simpleVector.h"
#include "modelCoeff.h"
#include "nameTools.h"
#include "indexTools.h"

// Helper class to check if the variable pattern in the model-term is
// compatible with modelMatrix needs.
class modelMatrixVariableChecker {
public:
   modelMatrixVariableChecker(parsedModelTerm & modeldescr) {
      OK = false;
      // also check if there are two variables ...
      if(modeldescr.variableTypes.size() != 2)
         throw generalRbayzError("modelMatrix requires two variables: a factor and a matrix.");
      bool acceptable0VarType = modeldescr.variableTypes[0]==1 || modeldescr.variableTypes[0]==2 ||
                              modeldescr.variableTypes[0]==4 || modeldescr.variableTypes[0]==5;
      if( acceptable0VarType && modeldescr.variableTypes[1]==6 )
         OK = true;
      else
         throw generalRbayzError("variable types in rr() model are not (convertable to) <factor>/<matrix>");

   }
   bool OK;
};

// Now modelMatrix can include a variableCheck member and initialize it in the constructor with
// the model-description. If that passes without throwing an error, the pattern is OK.
// If needed the status of the variable check can always be accessed as variableCheck.OK, but if
// no error was thrown it should have value True.

class modelMatrix : public modelCoeff {

/* working here ...
I'm trying to change F and M to member variables, which are then nicest initialized in the
initialization list (if not I would have to add 'empty' constructors and add methods to set them later).
To allow use in the initialization list, the checks done in the constructor may need to move to the
main bayz() function before calling these object constructors.... (is now done ...)
I'm seeing a disadvantage of moving it to main bayz(), this check is generic for Rreg, Lasso, Mixt
and all other cases that could work on a matrix of covariates. When these are split on the model allocation
tree, the check may need to be repeated when done in the main function. When staying here and all these
cases you modelMatrix as framework, it can be done only once here.
This then may need some checking function / class, that can also go in the initialization list....
*/
public:
   
   modelMatrix(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod)
   {
      // For now only allowing a matrix input where there is an index variable (model
      // made with id/matrix). It could be extended to allow for no id, so that matrix needs to
      // be aliged 1:1 with data records, then the 'id' is bascially a 1:1 link.
      F = new dataFactor(modeldescr.variableObjects[0], modeldescr.variableNames[0]);
      M = new dataMatrix(modeldescr.variableObjects[1], modeldescr.variableNames[1]);
      par = new parVector(modeldescr, 0.0l, M->colnames);
//      weights.initWith(M->ncol,1.0l);  // I think weights is not used (but using varmodel->weights)
      builObsIndex(obsIndex,F,M);
   }
   
   ~modelMatrix() {
      delete M;
      delete F;
      delete par;
   }
   
   // methods for single-site updates: data (de)corrections for a single covariate
   // column, and LHS and RHS statistics for a single covariate column.
   // Updated to skip (de)corrections when regression coeff is zero - this could have big
   // impact using the code for mixture models with many zero regcoeff.
   void resid_correct(size_t col) {
      if(par->val[col]==0.0l) return;
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->nelem; obs++) {
         resid[obs] -= par->val[col] * colptr[obsIndex[obs]];
      }
   }

   void resid_decorrect(size_t col) {
      if(par->val[col]==0.0l) return;
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->nelem; obs++)
         resid[obs] += par->val[col] * colptr[obsIndex[obs]];
   }

   // de+correct residuals and fit for a 'beta update': a change in beta.
   // The difference old-beta minus new-beta is passed as 'beta_diff'.
   void resid_fit_betaUpdate(double beta_diff, size_t col) {
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->nelem; obs++) {
         resid[obs] += beta_diff * colptr[obsIndex[obs]];
         fit[obs] -= beta_diff * colptr[obsIndex[obs]];
      }
   }

   // adjust residuals and fit for a change in scale, this works on total fit
   // (no argument for the column of the covar to use).
   // Notes:
   // 1) This could 'move up' to modelCoeff as it could be a generic function to
   //    update variances (obsviously only for random effects, but probably useful
   //    for mixture, lasso and ridge models).
   // 2) Now I'm doubting of the need to adjust fit-vector, it's like a known covariate
   //    here where we update a regression on that covariate, and that should not modify
   //    the covariate itself. But it can be this code has been 'harmless' so far because
   //    this fit is not used in other places and will be rebuilt in modelMatrix implementations.
   void resid_fit_scaleUpdate(double scale_old, double scale_new) {
      double resid_adjust = (1.0l - scale_new/scale_old);
      double fit_adjust = scale_new/scale_old;
      for (size_t obs=0; obs < F->nelem; obs++) {
         resid[obs] += fit[obs]*resid_adjust;
         fit[obs] = fit[obs]*fit_adjust;
      }
   }

   // [ToDo] I thought there could be some efficiency gain making a method that combined resis_decorrect()
   // and collect_lhs_rhs() because it is the same loops and some computations are re-used, but it looks more
   // difficult to take advantage of skipping zero regcoeff (which also partly happens in sample()) ....
   void collect_lhs_rhs(double & lhs, double & rhs, size_t col) {
      size_t matrixrow;
      lhs = 0.0l; rhs = 0.0l;
      double * colptr = M->data[col];
      double temp1;
      for (size_t obs=0; obs < F->nelem; obs++) {
         matrixrow = obsIndex[obs];
         temp1 = colptr[matrixrow] * residPrec[obs];
         rhs += temp1 * resid[obs];
         lhs += temp1 * colptr[matrixrow];
      }
   }
   
   // Note sure what I was making this one for, it is not (yet) used anywhere.
   void collect_sse(double & sse) {
      sse=0.0l;
      for (size_t obs=0; obs < F->nelem; obs++)
         sse += resid[obs]*resid[obs]*residPrec[obs];
   }

   // [ToDo] OBS: this fillFit is not right if the variance of random effect
   // is modelled with a scale. In that case fit = scale * par[] * covar[]
   void fillFit() {
      double * colptr;
      for (size_t obs=0; obs < F->nelem; obs++) fit[obs] = 0.0l;
      for(size_t k=0; k < M->ncol; k++) {
         colptr = M->data[k];
         for (size_t obs=0; obs < F->nelem; obs++)
            fit[obs] += par->val[k] * colptr[obsIndex[obs]];
      }
   }

   // Here no sample() yet, modelMatrix remains virtual. The derived classes implement sample()
   // by combining update_regressions() with update of hyper-paramters for that derived class.

   dataMatrix *M;
   dataFactor *F;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<size_t> obsIndex;

};

#endif /* modelMatrix_h */
