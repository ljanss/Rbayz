//
// BayzR -- dataAlgorithm.h
//
// Base class for dataAlgorithm classes - classes that define
// computational algorithms (MCMC related) on the different data classes.
// This mainly defines common methods; the data and helper vectors are very
// particular for the different algorithms. 
//
// Created by Luc Janss on 09/10/2026.

#ifndef dataAlgorithm_h
#define dataAlgorithm_h

#include "simpleVector.h"

class dataAlgorithm {
   
public:

   // There is no constructor defined here, there is no data in the base class, and the
   // derived classes have quite different data members and initialization requirements.
   // Generally, the constructors of the derived classes take data containers as arguments,
   // but the type and number differs.

   // But a virtual destructor can be defined.
   virtual ~dataAlgorithm() = default;

   // Two pure virtual versions of fillFit, one with an extra scale parameter, so all
   // derived classes must implement both versions! This kind of makes sense because some
   // coefficient models may need one, and other may need the other, even for the same data class.
   virtual void fillFit(simpleDblVec& fit, simpleDblVector& beta) = 0;
   virtual void fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale) = 0;

   // Collect LHS and RHS statistics; to make a common interface this is not defined to
   // use lhs and rhs vectors (in a few cases such as intercept and simple regression it is
   // a single element only).
   virtual void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs) = 0;

   // I have not yet considered putting residual correction/decorrection methods in the base
   // class here, mainly because some work by column (when working a dataMatrix structures), and
   // unlike fillFit is makes less sense that all derived classes should need to implement two
   // versions, or that the common interface has a 'column' argument that is not used in many cases.
};

#endif /* dataAlgorithm_h */