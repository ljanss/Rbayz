// modelCoeff.h
// Parent of all coefficient models (all fix/ran/reg fitting objects).
// Created by Luc Janss on 08/03/2021.

#ifndef modelCoeff_h
#define modelCoeff_h

#include <Rcpp.h>
#include "modelResp.h"
#include "parsedModelTerm.h"
#include <unistd.h>

class modelCoeff : public modelBase {
   
public:

   modelCoeff(parsedModelTerm & modeldescr, modelResp * rmod);

   ~modelCoeff() {
   }

   // The modelCoeff class has some methods that can be used to fit variance as
   // a scaling factor. These methods do not rely on the underlying data and only
   // use residuals and fitted values. These methods therefore are best placed here
   // because modelCoeff classes are the classes that have access to both resid and fit.
   void getFitScaleStats(double & lhs, double & rhs);
   void resid_fit_scaleUpdate(double scale_old, double scale_new);

   // Note sure what I was making this one for, it is not (yet) used anywhere, but
   // it also is a general method that could belong to modelCoeff, or else modelResp.
   void collect_sse(double & sse);

   modelResp* respModel;
   // Use of resid, residPrec and Nresid is for convenience so that all modelCoeff objects
   // have easy access to these vectors and their size. It only needs to be set once
   // here in modelCoeff ctor to be uniformly available in all derived classes.
   // It is implemented as a reference to the 'Values' member of the residuals (which are
   // actually parVector type), and a reference to the 'weights' of the variance model.
   simpleDblVector& resid;
   simpleDblVector& residPrec;
   size_t Nresid=0;
   simpleDblVector fit;

};

#endif /* modelCoeff_h */
