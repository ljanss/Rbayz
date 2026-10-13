//
// BayzR -- dataAlgoIntercept.h
//
// Algorithms to collect statistics for updating intercept (mean).
//
// Created by Luc Janss on 13/10/2026.

#ifndef dataAlgoIntercept_h
#define dataAlgoIntercept_h

#include "simpleVector.h"
#include "dataAlgorithm.h"

class dataAlgoIntercept : public dataAlgorithm {
   
public:

   dataAlgoIntercept() = default;

   ~dataAlgoIntercept() = default;

   void fillFit(simpleDblVec& fit, simpleDblVector& beta);

   void fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale);

   void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, double& lhs, double& rhs);

};

#endif /* dataAlgoIntercept_h */