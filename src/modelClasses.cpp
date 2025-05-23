//  BayzR --- modelClasses.cpp
//

#include <Rcpp.h>
#include "model_rn_cor.h"
#include "indexTools.h"
#include "optionsInfo.h"

/* *********************** model_rn_cor classes ***************************/

// ---- model_ran_cor_k0

// rn_cor_k0 can run 1 kernel or multiple kennels with "mergedKernel" appraoach - however, at the
// moment no option in the interface to toggle merging kernels or not.
model_rn_cor_k0::model_rn_cor_k0(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod)
{
   // for rn_cor_k0 all variance objects must be kernels (with an attached RObject)
   std::vector<varianceSpec> varianceList = modeldescr.allOptions.Vlist();
   for(size_t i=0; i<varianceList.size(); i++) {
      if (! varianceList[i].iskernel ) {   // iskernel true guarantees there is an RObject
         throw(generalRbayzError("Attempting to run rn_cor_k0 with parameterised kernels"));
      }
   }
   // Get the first kernel and then add (making kronecker products) with second etc., if available
   kernelList.push_back(new kernelMatrix(varianceList[0]));
   if (varianceList.size()==2) {  // combine with a second kernel if present
      kernelMatrix* K2 = new kernelMatrix(varianceList[1]);
      kernelList[0]->addKernel(K2);
      delete K2;
   }
   if (varianceList.size()>2) {  // need to think if I can keep combining kernels with addKernel()
      throw(generalRbayzError("Not yet ready to combine more than 2 kernels for interaction"));
   }
   // Here add a vector regcoeff (size K->ncol) to hold the regresssion on eigenvectors.
   // It is a parVector class so that the variance object can accept and work on it.
   regcoeff = new parVector(modeldescr, 0.0l, kernelList[0]->colnames);
   // obsIndex makes new level codes matching F->labels from every row in data to K->labels, it could
   // in principle replace the F->data and no need for the obsIndex vector.
   builObsIndex(obsIndex,F,kernelList[0]);
   // [ToDo] create the variance object - may need to move out as in ranfi
   varmodel = new diagVarStr(modeldescr, this->regcoeff, kernelList[0]->weights);
}


model_rn_cor_k0::~model_rn_cor_k0() {
   for(size_t i=0; i< kernelList.size(); i++)
      delete kernelList[i];
   delete regcoeff;
   delete varmodel;
}


void model_rn_cor_k0::sample() {
   // Update regressions on the eigenvectors
   double lhsl, rhsl;
   size_t matrixrow;
   double* colptr;
   kernelMatrix* K = kernelList[0];
   for (size_t obs=0; obs < F->nelem; obs++)
      fit.data[obs] = 0.0l;
   for(size_t col=0; col < K->ncol; col++) {
      colptr = K->data[col];
      // residual de-correction for this evec column
      for (size_t obs=0; obs < F->nelem; obs++)
         resid[obs] += regcoeff->val[col] * colptr[obsIndex[obs]];
      // Make the lhs and rhs and update this column regression
      lhsl = 0.0l; rhsl=0.0l;
      for (size_t obs=0; obs < F->nelem; obs++) {
         matrixrow = obsIndex[obs];
         rhsl += colptr[matrixrow] * residPrec[obs] * (resid[obs]-fit.data[obs]);
         lhsl += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
      }

      lhsl += varmodel->weights[col];
      regcoeff->val[col] = R::rnorm( (rhsl/lhsl), sqrt(1.0/lhsl));
      // residual correction for this column with updated regression
      for (size_t obs=0; obs < F->nelem; obs++)
         fit.data[obs] += regcoeff->val[col] * colptr[obsIndex[obs]];
   }
   for (size_t obs=0; obs < F->nelem; obs++)
      resid[obs] -= fit.data[obs];
}

void model_rn_cor_k0::sampleHpars() {
   varmodel->sample();
}

void model_rn_cor_k0::restart() {
   varmodel->restart();
}

// fillFit() here defines an empty version - making fit is already done in sample()
void model_rn_cor_k0::fillFit() { }

// prepForOutput puts the transform to random effects in the par-vector
// [!] this now only for 1, or 1 merged, kernel.
// Also: is this not the same as the fitted value??
void model_rn_cor_k0::prepForOutput() {
   kernelMatrix* K = kernelList[0];                  // kernelList is a std::vector<kernelMatrix*>
   for(size_t row=0; row< K->nrow; row++) {
      par->val[row]=0.0l;
      for(size_t col=0; col<K->ncol; col++) {
         par->val[row] += K->data[col][row] * regcoeff->val[col];
      }
   }
};

// ----- model_ran_cor_k1

/* start on making the "non-merged" approach where kernels remain individually in a list
model_rn_cor_k1::model_rn_cor_k1(parsedModelTerm & modeldescr, modelResp * rmod)
           : modelFactor(modeldescr, rmod), regcoeff(), fitval(), gprior(modeldescr.allOptions["prior"]) {
   // For the moment all variance objects must be kernels
   for(size_t i=0; i<modeldescr.varObject.size(); i++) {
      if (modeldescr.varObject[i]==R_NilValue) {
         throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
      }
   }
   for(size_t i=0; i<modeldescr.varObject.size(); i++) {
      kernelList.push_back(new kernelMatrix(modeldescr.varObject[i], modeldescr.varName[i]));
   }
}
*/

