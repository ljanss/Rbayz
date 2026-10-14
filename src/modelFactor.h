//
//  rbayz -- modelFactor.h
//  Computational methods to work on one factor, but that can be coming from a model-term with interaction
//  that is recoded in dataFactor to be stored again as one factor.
//  All modelFixf and modelRanfi use this set-up and the methods in this class for sampling / udpating.
//  From the correlated random models, Ranfc1 uses the data set-up form this class, but not the methods.
//  This is still not a concrete class! (Has no sample() implemented).
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelFactor_h
#define modelFactor_h

#include <Rcpp.h>
#include <cmath>
#include "dataAlgorithm.h"
#include "modelCoeff.h"
#include "dataFactor.h"
#include "optionsInfo.h"
#include "simpleVector.h"
//#include <unistd.h>

class modelFactor : public modelCoeff {

public:
   
   modelFactor(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod),
           F(modeldescr.variableObjects, modeldescr.variableNames),
           comp(F),
           lhs(F.labels.nelem,0),
           rhs(F.labels.nelem,0)
   {
      par = new parVector(modeldescr, 0.0l, F.labels);
   }

   // constructor with a variance list used by (some) random effect models; it is used to
   // take labels for the factors if there are kernels in the variance list.
   modelFactor(parsedModelTerm & modeldescr, modelResp * rmod, std::vector<varianceSpec> varlist)
         : modelCoeff(modeldescr, rmod),
           F(modeldescr.variableObjects, modeldescr.variableNames, varlist),
           comp(F),
           lhs(F.labels.nelem,0),
           rhs(F.labels.nelem,0)
   {
      par = new parVector(modeldescr, 0.0l, F.labels);
   }
   
   ~modelFactor() {
      delete par;
   }

   dataFactor F;
   dataAlgoFactor comp;
   simpleDblVector lhs, rhs;

};

#endif /* modelFactor_h */
