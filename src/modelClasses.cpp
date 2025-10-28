//  BayzR --- modelClasses.cpp
//

#include <Rcpp.h>
#include "Rbayz.h"
#include "model_rn_cor.h"
#include "indexTools.h"
#include "optionsInfo.h"

/* *********************** model_rn_cor classes ***************************/

// ---- model_ran_cor_k0

// rn_cor_k0 runs models with 1 kernel or multiple kernels.
// With "mergedKernel", multiple kernels are merged into one and the code runs on one kernel,
// otherwise it will apply the 'on the fly' construction of the kernel kronecker products. 
model_rn_cor_k0::model_rn_cor_k0(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod, modeldescr.allOptions.Vlist(), modeldescr.allOptions["mergeKernels"].valbool), 
        alpha2levels(), regcoeff(nullptr), varmodel(nullptr)
{

   std::vector<varianceSpec> varianceList = modeldescr.allOptions.Vlist();
   if(varianceList.size() != F->Nvar) {
      throw(generalRbayzError("The number of interaction variables in [" + modeldescr.shortModelTerm + "] does not match the number of variance terms"));
   }

   // for rn_cor_k0 all variance objects must be kernels (with an attached RObject)
   for(size_t i=0; i<varianceList.size(); i++) {
      if (! varianceList[i].iskernel ) {   // iskernel true guarantees there is an RObject
         throw(generalRbayzError("Attempting to run rn_cor_k0 with parameterised kernels"));
      }
   }

   // Check model term vdimp option; this is used to reset the default dimp=90 to select evecs in each kernel.
   // Note: I was consdidering to also allow a vdim, but that's not yet implemented, and the current kernelMatrix
   // constructor can only be tuned on 'dimp' but not on 'dim'.
   if(modeldescr.allOptions["vdimp"].isgiven) {
      var_retain = modeldescr.allOptions["vdimp"].valnumb[0];
      if(var_retain < 10 || var_retain > 100.0)
         throw(generalRbayzError("vdimp option should be between 10 and 100"));
      var_retain = pow(var_retain/100.0, (1.0/double(varianceList.size())))*100.0; // convert to per-kernel pct
   }
   else {
      var_retain = 90;
   }

   // Store all kernels in the kernelList vector; this includes doing the eigen-decomposition
   // and selecting number of eigenvectors to retain in each. 
   for(size_t i=0; i<varianceList.size(); i++) {
      kernelList.push_back(new kernelMatrix(varianceList[i], var_retain));
   }

   // Check sizes if kernels would be combined by kronecker product; the merged_ncol will also
   // be the size of the alpha vector to be estimated.
   size_t merged_nrow=1, merged_ncol=1;
   for(size_t i=0; i< kernelList.size(); i++) {
      merged_nrow *= kernelList[i]->nrow;
      merged_ncol *= kernelList[i]->ncol;
   }

   // set maxmen: max memory to use for large matrices. If there is no setting in the model-term, set to 4GB.
   // A maxmem option from the model-term should be in GB, so multiply by 1e9 to get bytes.
   size_t maxmem;
   if (modeldescr.allOptions["maxmem"].isgiven) {
      maxmem = (size_t) (1e9 * modeldescr.allOptions["maxmem"].valnumb[0] );
   }
   else
      maxmem = 4e9;  // default maxmem = 500 million doubles = 4 GB

   // Check to merge kernels; for now only do this when indicated by the user.
   if (modeldescr.allOptions["mergeKernels"].isgiven && modeldescr.allOptions["mergeKernels"].valbool) {
      if (varianceList.size() < 2) {
         Rbayz::Messages.push_back("Warning: mergeKernels option for <" + modeldescr.shortModelTerm + "> is ignored, there is only one kernel");
      } else {
         size_t mem_needed = merged_nrow * merged_ncol * 8; // in bytes, double = 8 bytes
         if( mem_needed > maxmem ) {
            throw(generalRbayzError("Merging kernels needs " + std::to_string((size_t)(mem_needed/1e9)+1) + "GB, increase maxmem or do not merge kernels"));
         }
         // if here, memory is OK, so merge all kernels into the first one in the list
         for(size_t i=1; i< kernelList.size(); i++) {
            kernelList[0]->addKernel(kernelList[i]);
            delete kernelList[i];
            kernelList.erease(kernelList.begin()+i);
         }
      }
      // Now merged_ncol across all must have become dimension of the first kernel
      if( varianceList.size() != 1 || merged_ncol != kernelList[0]->ncol ) {
         throw(generalRbayzError("Something went wrong merging kernels, please consult the developers"));
      }
   }

   // If multiple kernels are not merged, we need some mapping of combinations of evecs in the kernels
   // to elements in the regcoeff vector. This is done in the alpha2levels matrix.
   // Note: the factors in the interaction model is put on rows, the alpha's on columns.
   if (kernelList.size() > 1) {
      if( merged_ncol > 100000 ) {
         // [ToDo] If messages could have an option that displays messages directly on the screen, this could be a useful one to show up immediately.
         Rbayz::Messages.push_back("Warning: the number of regressions modeled in <" + modeldescr.shortModelTerm + "> is large (" + std::to_string(merged_ncol) + ")";
      }
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

   // Set-up / initialize 'regcoeff' (alpha) regression vector; for multiple kernels need to prepare
   // labels; for one kernel it can use the column names from kernel[0].
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
   else { // only one kernel
      regcoeff = new parVector(modeldescr, 0.0l, kernelList[0]->colnames);
   }
   regcoeff->Name = regcoeff->Name + ".alpha";
   if(modeldescr.allOptions["alpha_est"].isgiven || modeldescr.allOptions["alpha_save"].isgiven) {
      Rbayz::parList.push_back(&regcoeff);
   }
   if (modeldescr.allOptions["alpha_save"].isgiven) {
      regcoeff->saveSamples=true;
   }

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

   // [ToDo]: sample() needs to be updated with separate code to run on one or on multiple kernels.
   // For one kernel, a variable like alpha2levels is not initialised and should not be used!

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

