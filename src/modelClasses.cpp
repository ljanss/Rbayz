//  BayzR --- modelClasses.cpp
//

#include <Rcpp.h>
#include "Rbayz.h"
#include "model_rn_cor.h"
#include "indexTools.h"
#include "optionsInfo.h"

/* ****************** modelRanfc classes ***************************/

// --------------------- modelRanfc1 ------------------------------

/* For the Ranfc1 version, interactions and kernels should be recoded / merged so code can work
   as if handling one factor and one kernel (including the case if there was really only one
   factor with one kernel). If it originates from a model-term with interaction and multiple kernels,
   it should have the 'mergeKernels' options.
*/
// rn_cor_k0 runs models with 1 kernel or multiple kernels.
// With "mergedKernel", multiple kernels are merged into one and the code runs on one kernel,
// otherwise it will apply the 'on the fly' construction of the kernel kronecker products. 
modelRanfc1::modelRanfc1(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod, modeldescr.allOptions.Vlist()), 
        regcoeff(nullptr), varmodel(nullptr)
{

   std::vector<varianceSpec> varianceList = modeldescr.allOptions.Vlist();
   if(varianceList.size() != F->Nvar) {
      throw(generalRbayzError("The number of interaction variables in [" + modeldescr.shortModelTerm + "] does not match the number of variance terms"));
   }

   // for Ranfc1 all variance objects must be kernels (with an attached RObject)
   for(size_t i=0; i<varianceList.size(); i++) {
      if (! varianceList[i].iskernel ) {   // iskernel true guarantees there is an RObject
         throw(generalRbayzError("Error: running Ranfc1 with parameterised kernels - pls report to developers"));
      }
   }

   // If there are multiple kernels, the model should have 'mergeKernels' option.
   // The default is now NOT merging. To merge, the option should be given, and should be true.
   // Note: the right triage should have been done in rbayz main, so ending up here without
   // mergeKernels options is a programming bug.
   if (varianceList.size() > 1) {
      if (!(modeldescr.allOptions["mergeKernels"].isgiven && modeldescr.allOptions["mergeKernels"].valbool)) {
         throw generalRbayzError("Error: running Ranfc1 without merging kernels - pls report to developers");
      }
   }

   // Check model term vdimp option; this is used to reset the default dimp=90 to select evecs in each kernel.
   // Note: I was consdidering to also allow a vdim, but that's not yet implemented, and the current kernelMatrix
   // constructor can only be tuned on 'dimp' but not on 'dim'.
   double var_retain;
   if(modeldescr.allOptions["vdimp"].isgiven) {
      var_retain = modeldescr.allOptions["vdimp"].valnumb[0];
      if(var_retain < 10 || var_retain > 100.0)
         throw(generalRbayzError("vdimp option should be between 10 and 100"));
      if (varianceList.size() > 1)
         var_retain = pow(var_retain/100.0, (1.0/double(varianceList.size())))*100.0; // convert to per-kernel pct
   }
   else {
      var_retain = 90;
   }

   // set maxmen: max memory to use for large matrices. If there is no setting in the model-term, set to 4GB.
   // A maxmem option from the model-term should be in GB, so multiply by 1e9 to get bytes.
   size_t maxmem;
   if (modeldescr.allOptions["maxmem"].isgiven) {
      maxmem = (size_t) (1e9 * modeldescr.allOptions["maxmem"].valnumb[0] );
   }
   else
      maxmem = 4e9;  // default maxmem = 500 million doubles = 4 GB

   if (varianceList.size() == 1) { // one kernel, allocate directly to the 'kernel' member variable
      K = new kernelMatrix(varianceList[0], var_retain);
   }
   else { // multiple kernels, need merging
      // 1. Store all kernels in the kernelList vector; this includes doing the eigen-decomposition
      // and selecting number of eigenvectors to retain in each. kernelList is local for temporary use.
      std::vector<kernelMatrix*> kernelList;
      for(size_t i=0; i<varianceList.size(); i++) {
         kernelList.push_back(new kernelMatrix(varianceList[i], var_retain));
      }
      // 2. Check sizes if kernels would be combined by kronecker product; the merged_ncol will also
      // be the size of the alpha vector.
      size_t merged_nrow=1, merged_ncol=1;
      for(size_t i=0; i< kernelList.size(); i++) {
         merged_nrow *= kernelList[i]->nrow;
         merged_ncol *= kernelList[i]->ncol;
      }
      if( merged_ncol > 100000 ) {
         // [ToDo] If messages could have an option that displays messages directly on the screen,
         // this could be a useful one to show up immediately.
         Rbayz::Messages.push_back("Warning: the number of regressions modeled in <" + modeldescr.shortModelTerm +
             "> is large (" + std::to_string(merged_ncol) + ")";
      }
      size_t mem_needed = merged_nrow * merged_ncol * 8; // in bytes, double = 8 bytes
      if( mem_needed > maxmem ) {
         throw(generalRbayzError("Merging kernels needs " + std::to_string((size_t)(mem_needed/1e9)+1) + "GB, increase maxmem or do not merge kernels"));
      }
      // if here, memory is OK, so merge all kernels into the first one in the list
      for(size_t i=1; i< kernelList.size(); i++) {
         kernelList[0]->addKernel(kernelList[i]);
      }
      // Now merged_ncol across all must have become dimension of the first kernel
      if( merged_ncol != kernelList[0]->ncol ) {
         throw(generalRbayzError("Something went wrong merging kernels, please consult the developers"));
      }
      // if all OK and done, copy/swap kernelList[0] to the 'K' member variables.
      // Note: swap is from simpleMatrix parent class, so kernelList[0] in the swap call will decast the
      // kernelMatrix to a simpleMatrix, and only work on copying the 'bare' matrix data.
      // [ToDo]: it could be nice if labeledMatrix and kernelMatrix also have a swap, but they should
      // work with the hierarchy. Then each class swap() adds copying its own added member variables ... hmmm
      // it looks like such a thing can't work because a derived class swap overrides the parent swap. 
      // But ... I search on making copy contructors in a hierarchy, which should be possible. 
      K.swap(kernelList[0]);
      K.rownames = kernelList[0].rownames;    // members from labeledMatrix, these are std::vector and
      K.colnames = kernelList[0].colnames;    // will copy nicely ..
      K.weights.swap(kerneList[0].weights);   // weights from kernelMatrix is a simpleDblVector, which also has swap()
      K.sumEvalues = kernelList[0].sumEvalues;
      // clean up 'new' allocations - the kernelList vector will clean up itself.
      for(size_t i=0; i< kernelList.size(); i++)
         delete kernelList[i];
   }

   // Set-up / initialize 'regcoeff' (alpha) regression vector
   regcoeff = new parVector(modeldescr, 0.0l, K->colnames);
   regcoeff->Name = regcoeff->Name + ".alpha";
   if(modeldescr.allOptions["alpha_est"].isgiven || modeldescr.allOptions["alpha_save"].isgiven) {
      Rbayz::parList.push_back(&regcoeff);
   }
   if (modeldescr.allOptions["alpha_save"].isgiven) {
      regcoeff->saveSamples=true;
   }

   // obsIndex makes new level codes matching F->labels from every row in data to K->labels, it could
   // in principle replace the F->data and no need for the obsIndex vector.
   // I'm leaving it for a moment to see in debugging if obsIndex is indeed the same as F->data,
   // in the sample code it is already removed.
   builObsIndex(obsIndex,F,kernelList[0]);

   // [ToDo] create the variance object - may need to move out as in ranfi when allowing for
   // different variance structures. But this is the variance structure for the alpha coefficients,
   // and there is no interface yet to allow different structures here...
   varmodel = new diagVarStr(modeldescr, this->regcoeff, kernelList[0]->weights);
}

modelRanfc1::~modelRanfc1() {
   for(size_t i=0; i< kernelList.size(); i++)
      delete kernelList[i];
   delete regcoeff;
   delete varmodel;
}

void modelRanfc1::sample() {
   // Update regressions on the eigenvectors
   double lhsl, rhsl;
   size_t matrixrow;
   double* colptr;
   for (size_t obs=0; obs < F->nelem; obs++)
      fit.data[obs] = 0.0l;
   for(size_t col=0; col < K->ncol; col++) {
      colptr = K->data[col];
      // residual de-correction for this evec column
      for (size_t obs=0; obs < F->nelem; obs++)
         resid[obs] += regcoeff->val[col] * colptr[F->data[obs]];
      // Make the lhs and rhs and update this column regression
      lhsl = 0.0l; rhsl=0.0l;
      for (size_t obs=0; obs < F->nelem; obs++) {
         matrixrow = F->data[obs];
         rhsl += colptr[matrixrow] * residPrec[obs] * (resid[obs]-fit.data[obs]);
         lhsl += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
      }

      lhsl += varmodel->weights[col];
      regcoeff->val[col] = R::rnorm( (rhsl/lhsl), sqrt(1.0/lhsl));
      // residual correction for this column with updated regression
      for (size_t obs=0; obs < F->nelem; obs++)
         fit.data[obs] += regcoeff->val[col] * colptr[F->data[obs]];
   }
   for (size_t obs=0; obs < F->nelem; obs++)
      resid[obs] -= fit.data[obs];
}

void modelRanfc1::sampleHpars() {
   varmodel->sample();
}

void modelRanfc1::restart() {
   varmodel->restart();
}

// fillFit() here defines an empty version - making fit is already done in sample()
void modelRanfc1::fillFit() { }

// prepForOutput puts the transform to random effects in the par-vector
// [!] this now only for 1, or 1 merged, kernel.
// Also: is this not the same as the fitted value??
void modelRanfc1::prepForOutput() {
   for(size_t row=0; row< K->nrow; row++) {
      par->val[row]=0.0l;
      for(size_t col=0; col<K->ncol; col++) {
         par->val[row] += K->data[col][row] * regcoeff->val[col];
      }
   }
};

// ------------------------- modelRanfck ----------------------------

/* Ranfck is the class that handles multiple kernels, and only kernels (no US or other included).
   Cases with a single kernel (including 'merged' kernels) should have been sent
   by rbayz main to the Ranfc1 class.
*/

modelRanfck::modelRanfck(parsedModelTerm & modeldescr, modelResp * rmod)
           : modelCoeff(modeldescr, rmod), regcoeff(nullptr), gprior(modeldescr.allOptions["prior"]) {
   // For the moment all variance objects must be kernels
   for(size_t i=0; i<modeldescr.varObject.size(); i++) {
      if (modeldescr.varObject[i]==R_NilValue) {
         throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
      }
   }
   for(size_t i=0; i<modeldescr.varObject.size(); i++) {
      kernelList.push_back(new kernelMatrix(modeldescr.varObject[i], modeldescr.varName[i]));
   }

   // If multiple kernels are not merged, we need some mapping of combinations of evecs in the kernels
   // to elements in the regcoeff vector. This is done in the alpha2levels matrix.
   // Note: the factors in the interaction model is put on rows, the alpha's on columns.
   if (kernelList.size() > 1) {
      alpha2levels = simpleIntMatrix(kernelList.size(), merged_ncol);
      size_t prev_levs, this_levs, next_levs;
      for(size_t i=0; i< kernelList.size(); i++) {
         this_levs = F->factorList[i]->nlevels;
         prev_levs = 1; next_levs = 1;
         for(size_t j=0; j<i; j++)
            prev_levs *= F->factorList[j]->nlevels;
         for(size_t j=i+1; j< kernelList.size(); j++)
            next_levs *= F->factorList[j]->nlevels;
         for(size_t k1=0; k1< prev_levs; k1++) {
            for(size_t k2=0; k2< this_levs; k2++) {
               for(size_t k3=0; k3< next_levs; k3++) {
                  size_t col = k1*this_levs*next_levs + k2*next_levs + k3;
                  alpha2levels.data[i][col] = k2;
               }
            }
         }
      }
   }

      if (kernelList.size() > 1) {
      std::vector<std::string> temp_labels(merged_ncol);
      for(size_t col=0; col< merged_ncol; col++) {
         std::string nm = kernelList[0]->colnames[ alpha2levels.data[0][col] ];
         for(size_t row=1; row< kernelList.size(); row++) {
            nm += "." + kernelList[row]->colnames[ alpha2levels.data[row][col] ];
         }
         temp_labels[col] = nm;
      }
      regcoeff = new parVector(modeldescr, 0.0l, temp_labels);
   }

}

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

