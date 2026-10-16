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

   virtual void fillFit() = 0;

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
   // The following is for convenience so that all modelCoeff objects have direct
   // pointers to residuals and residual variance, and it only needs to be set once
   // here in modelCoeff cstr'or.
   // [ToDo] It is a bit old-style C that these resid and residPrec are raw pointers and
   // then used as arrays, it would be nicer if this was e.g. a reference to simpleDblVector.
   simpleDblVector& resid;
   double *residPrec=NULL;
   size_t Nresid=0;
   simpleDblVector fit;

};

#endif /* modelCoeff_h */
