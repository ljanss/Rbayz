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

   // Note: methods to collect LHS and RHS statistics have different interfaces in different
   // dataAlgorithm classes and are therefore not part of the base class. 
   // Note: I'm also no longer including residual correct/decorrect methods because the code
   // gets updated to do this through adding/subtracting the fit-vectors, which can be done in one statement. 
//   virtual void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs) { }
//   virtual void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, double& lhs, double& rhs, size_t column) { }

};

#endif /* dataAlgorithm_h */