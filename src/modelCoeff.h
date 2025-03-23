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

   modelCoeff(parsedModelTerm & modeldescr, modelResp * rmod) : modelBase(), fit() {
      // updating bayzR to handle hierarchical models can imply some (large?) changes here,
      // because then the "response model" can be another coefficient model.
      // modelResp and modelCoeff are in different branches of the class hierarchy, and cannot
      // be cast from one to the other.
      // Maybe it needs a special interface class to help presenting a modelCoeff object
      // "as if" it is a modelResp class for the hierarchical model?
      respModel=rmod;
      resid = respModel->resid->val;
      if(respModel->varModel==0) Rcpp::Rcout << "!respModel->varModel pointer is zero!\n";
      residPrec = respModel->varModel->weights.data;
      Nresid = respModel->resid->nelem;
      fit.initWith(Nresid, 0.0l);
   }

   ~modelCoeff() {
   }

   virtual void fillFit() = 0;

   // make statistics to estimate scale of fitted values:
   // lhs=sum(fit^2), rhs=sum(fit*resid) with resid de-corrected for fit.
   void getFitScaleStats(double & lhs, double & rhs) {
      lhs=0.0l;
      rhs=0.0l;
      for(size_t i=0; i < Nresid; i++) {
         lhs += fit.data[i] * fit.data[i];
         rhs += fit.data[i] * (resid[i] + fit.data[i]);
      }
   }

   modelResp* respModel;
   // The following is for convenience so that all modelCoeff objects have direct
   // pointers to residuals and residual variance, and it only needs to be set once
   // here in modelCoeff cstr'or.
   double *resid=NULL, *residPrec=NULL;
   size_t Nresid=0;
   simpleDblVector fit;

};

#endif /* modelCoeff_h */
