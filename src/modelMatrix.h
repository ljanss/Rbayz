//
//  modelMatrix.h
//  rbayz
//
//  Defines the computational methods when design matrix is a (real) matrix:
//    -> has pointer to dataMatrix object
//    -> has pointer to dataFactor object
//  and defines residual de/correct and collect lhs/rhs methods working on this
//  kind of objects.
//  This is not yet a concrete class, derived classes differ mostly in the
//  constructors that define how different kinds of matrix data is prepared.
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

class modelMatrix : public modelCoeff {
   
public:
   
   modelMatrix(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod)
   {
      // For now only allowing a matrix input where there is an index variable (model
      // made with id/matrix). It could be extended to allow for no id, so that matrix needs to
      // be aliged 1:1 with data records, then the 'id' is bascially a 1:1 link.
      // 2: probably more variable-types could be accepted here, if dataFactor and dataMatrix
      // would know how to convert them (e.g. integer vector to factor, dataframe to matrix).
      if( ! (modeldescr.variableTypes[0]==1 && modeldescr.variableTypes[1]==6) )
         throw generalRbayzError("variable types in rr() model are not <factor>/<matrix>");
      F = new dataFactor(modeldescr.variableObjects[0]);
      M = new dataMatrix(modeldescr.variableObjects[1], modeldescr.variableNames[1]);
   // OBS M->colnames may not be filled, need to fix in labeledMatrix.
      par = new parVector(modeldescr.variableString, 0.0l, M->colnames);
      weights.initWith(M->ncol,1.0l);
      builObsIndex(obsIndex,F,M);
      lhs = 0.0l;
      rhs = 0.0l;
   }
   
   ~modelMatrix() {
      delete M;
      delete F;
   }
   
   void resid_correct(size_t col) {
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] -= par[col] * colptr[obsIndex[obs]];
   }

   void resid_decorrect(size_t col) {
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] += par[col] * colptr[obsIndex[obs]];
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t matrixrow;
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++) {
         matrixrow = obsIndex[obs];
         rhs += colptr[matrixrow] * residPrec[obs] * resid[obs];
         lhs += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
      }
   }

   void accumFit(simpleDblVector & fit) {
      // There is a fit vector declared in this object, but it is not filled.
      // If that fit would be available, addding it here to the total fit would be faster.
      double * colptr;
      for(size_t k=0; k < M->ncol; k++) {
         colptr = M->data[k];
         for (size_t obs=0; obs < F->data.nelem; obs++)
            fit[obs] += par[k] * colptr[obsIndex[obs]];
      }
   }

   // Here no sample() yet, modelMatrix remains virtual. The derived classes implement sample()
   // by combining update_regressions() with update of hyper-paramters for that derived class.

   /* old code from ranf_cor for comparison
   void resid_correct(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] -= par[col] * matrixdata(coldata(obs),col);
   }
   
   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] += par[col] * matrixdata(coldata(obs),col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t rowlevel;
      for (size_t obs=0; obs < coldata.size(); obs++) {
         rowlevel = coldata(obs);
         rhs += matrixdata(rowlevel,col) * residPrec[obs] * resid[obs];
         lhs += matrixdata(rowlevel,col) * matrixdata(rowlevel,col) * residPrec[obs];
      }
   }

   */
   dataMatrix *M;
   dataFactor *F;
   simpleDblVector weights;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   std::vector<size_t> obsIndex;

};

#endif /* modelMatrix_h */
