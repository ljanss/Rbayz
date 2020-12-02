//
//  BayzR --- modelRreg.hpp
//
//  Computational class to to model random regressions.
//  This used all methods from modelMatrix, only some parameter names need to be set.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRreg_h
#define modelRreg_h

#include <Rcpp.h>
#include "modelMatrix.h"
#include "dataMatrix.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelMatrix(modelTerm, d, e, resp)
   {
      // Note: modelMatrix checks and sets up factor and matrix, Rreg is an implementation
      // of modelMatrix that only needs to add modeling of variances
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelRreg() {
   }
   
   void sample() {
      update_regressions();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]/M->weights[k];
      hpar[0] = gprior.samplevar(ssq, M->ncol);
   }


};

#endif /* modelRreg */
