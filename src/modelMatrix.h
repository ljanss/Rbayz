// R/bayz
// modelMatrix.h
//
// Coefficient-model class that defines the skeleton of classes working with matrices with
// features (covariates). This is not yet a concrete class, but specifies that these coefficient
// model classses use two data containers (factor and matrix) and the dataAlgoMatrix computational
// methods.
// Derived concrete classes need to add prior specifics and sampling methods. Currently the
// derived class for rr() model-terms is implemented and working on mx(); potentially can also
// revive the (exponential) lasso in some form as a third sibling in this branch.
// Note that, although this class uses a dataFactor data container, it is not implemented as a
// derived class of modelFactor, mostly because an interaction of factors would mess up the data
// linking made here (a re-coded interaction does not exist as a single variable in the data).
// Here an interaction will not pass because the parser will detect this as three variables, which
// cannot be handled in the constructor.
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

// Helper class to check if the obtained variables are two variables (convertable to)
// type factor and matrix as needed in modelMatrix classes.
class modelMatrixVariableChecker {
public:
   modelMatrixVariableChecker(parsedModelTerm & modeldescr) {
      OK = false;
      if(modeldescr.variableTypes.size() != 2)
         throw generalRbayzError("Variables in " + modeldescr.shortModelTerm + 
               " seem not be of var1/var2 pattern.");
      bool acceptable0VarType = modeldescr.variableTypes[0]==1 || modeldescr.variableTypes[0]==2 ||
                              modeldescr.variableTypes[0]==4 || modeldescr.variableTypes[0]==5;
      if( acceptable0VarType && modeldescr.variableTypes[1]==6 )
         OK = true;
      else
         throw generalRbayzError("Variables in " + modeldescr.shortModelTerm + 
               " are not (convertable to) <factor>/<matrix>");
   }
   bool OK;
};

// Note: for now only allowing a matrix input where there is an index variable (model
// made with id/matrix). It could be extended to allow for no id, so that matrix needs to
// be aliged 1:1 with data records, then the 'id' is bascially a 1:1 link.
// F and M are now initialized in the initialization list

class modelMatrix : public modelCoeff {
public:
   modelMatrix(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod),
           variableCheck(modeldescr),
           F(modeldescr.variableObjects[0], modeldescr.variableNames[0]),
           M(modeldescr.variableObjects[1], modeldescr.variableNames[1]),
           comp(F, M)
   {
      par = new parVector(modeldescr, 0.0l, M.colnames);
   }
   
   ~modelMatrix() {
      delete par;
   }
   
   // The evec fit (but not yet rr-fit) works with resid_decorrect per column, then
   // accumulate fit per column, uses corrected residual as current resid minus so-far
   // accumualted fit, and in the end does resisual correction by subtractin the whole
   // accumated fit. This can work using an fillAccumFit(..., col) method that does the
   // accumulated fit per column (and does not need resid_correct).

   // Here no sample() yet, modelMatrix remains virtual. The derived classes implement sample()
   // by combining update_regressions() with update of hyper-paramters for that derived class.

   modelMatrixVariableChecker variableCheck;
   dataFactor F;
   dataMatrix M;
   dataAlgoMatrix comp;

};

#endif /* modelMatrix_h */
