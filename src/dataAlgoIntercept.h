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

   void resid_correct(simpleDblVec& resid, simpleDblVector& beta);

   void resid_decorrect(simpleDblVec& resid, simpleDblVector& beta);

   void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs);

};

#endif /* dataAlgoIntercept_h */