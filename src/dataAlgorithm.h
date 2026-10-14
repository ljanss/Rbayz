//
// BayzR -- dataAlgorithm.h
//
// Data algorithm classes: classes the define the computations on the interface of
// data containers and coefficient models, typically by supplying the statistics needed
// for sampling (and PCG updating too?). 
// 
// The base class is relatively empty but defines the requirement for fillFit
// methods to be implemented: all classes need to define a fillFit with and
// without scale; this is essential for correct working of residual de/corrections,
// which is done with the fit-vectors and don't account for scale directly.
//
// Concrete classes are needed for every combination of data containers that can
// arise in model terms. Thus a single covariate regression (one dataCovar container)
// needs different algorithms than a nested covariate regression (needs a dataFactor
// and a dataCovar container). But the nested covariate algorithm can be re-used across
// fixed nested regression and random slope models.
// The concrete classes usually implement a collectLhsRhs() method, but argument
// signatures differ, and some container-combination may need more or other methods.
//
// Created by Luc Janss on 09/10/2026.

#ifndef dataAlgorithm_h
#define dataAlgorithm_h

#include "simpleVector.h"
#include "dataFactor.h"
#include "dataCovar.h"
#include "dataMatrix.h"

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

//   virtual void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, simpleDblVec& lhs, simpleDblVec& rhs) { }
//   virtual void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, double& lhs, double& rhs, size_t column) { }

};

// Data algorithms for intercept fit. As an exception, this is not connected to a data
// container, the intercept is trivial and implied (bayz does not explicitly store a vector
// of 1's to work on).
class dataAlgoIntercept : public dataAlgorithm {
public:
   dataAlgoIntercept() = default;
   ~dataAlgoIntercept() = default;
   void fillFit(simpleDblVec& fit, simpleDblVector& beta);
   void fillFit(simpleDblVec& fit, simpleDblVector& beta, double scale);
   void collectLhsRhs(simpleDblVec& resid, simpleDblVec& residprec, double& lhs, double& rhs);
};

// Data algorithms for dataFactor. Note the dataFactor standardly re-codes interactions
// of factors to behave again as one factor, so this also covers multi-factor interactions.
// Note: there is also a "not collapsed" dataFactorNC, but I will see if that also needs its
// own algorithm class.
class dataAlgoFactor : public dataAlgorithm {
public:
   dataAlgoFactor(dataFactor& Fdata_) : Fdata(Fdata_) {}
   ~dataAlgoFactor() = default;
   void fillFit(simpleDblVector& fit, simpleDblVector& beta);
   void fillFit(simpleDblVector& fit, simpleDblVector& beta, double scale);
   void collectLhsRhs(simpleDblVector& resid, simpleDblVector& residprec, simpleDblVector& lhs, simpleDblVector& rhs);
private:
   dataFactor& Fdata;
};

// Data algorithms for a single covariate (one dataCovar data container).
// For a single covariate, beta must be a scalar and lhs and rhs too.
class dataAlgoCovar : public dataAlgorithm {
public:
   dataAlgoCovar(dataCovar& Cdata_) : Cdata(Cdata_) {}
   ~dataAlgoCovar() = default;
   void fillFit(simpleDblVector& fit, double& beta);
   void fillFit(simpleDblVector& fit, double& beta, double scale);
   void collectLhsRhs(simpleDblVector& resid, simpleDblVector& residprec, double& lhs, double& rhs);
private:
   dataCovar& Cdata;
};

// Can here now develop nested regression too, a case with two types of data containers.



#endif /* dataAlgorithm_h */