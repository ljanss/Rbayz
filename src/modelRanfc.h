//  BayzR --- modelRanfc.h
//  Computational classes to model RANdom Factors with Correlations (from using rn(..., V=K)).
//   - Ranfc1: for one kernel, or 'mergedKernel' and recoded interactions (so it computes again as if
//     one random effect with one kernel). This derives from modelFactor that can do the data management
//     of recoding interaction to one factor, but methods from modelFactor cannot be used and are redone.
//   - Ranfck: for multiple kernels that are not recoded into one factor and one merged kernel. This
//     one does not derive from modelFactor (but directly from modelCoeff) and uses a dataFactorNC to
//     manage the factor data as a vector of factors.
// Note: implementations are in modelClasses.cpp.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanfc_h
#define modelRanfc_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "kernelMatrix.h"
#include "priorClasses.h"
#include "parsedModelTerm.h"
#include "simpleMatrix.h"

class modelRanfc1 : public modelFactor {

public:

   modelRanfc1(parsedModelTerm & modeldescr, modelResp * rmod);
   ~modelRanfc1();
   kernelMatrix* K;
   parVector *regcoeff;
   std::vector<size_t> obsIndex;
   indepVarStr* varmodel;

};

class modelRanfck : public modelCoeff {

   modelRanfck(parsedModelTerm & modeldescr, modelResp * rmod);
   ~modelRanfck();
   std::vector<kernelMatrix*> kernelList;
   parVector *regcoeff;
   indepVarStr* varmodel;
   simpleIntMatrix alpha2levels;

}


#endif /* modelRanfc_h */
