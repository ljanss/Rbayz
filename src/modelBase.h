//  modelBase.h
//
//  Base class for model (computational) classes.
//  This defines common interface with par-vector and sample() vector, many
//  other details specific for response, explanatory variables or variances come
//  in derived classes.
//
//  Created by Luc Janss on 03/08/2018.

#ifndef modelBase_h
#define modelBase_h

#include <Rcpp.h>
#include <vector>
#include <stdio.h>
#include "Rbayz.h"
#include "parVector.h"
//#include <unistd.h>

class modelBase {
   
public:

   modelBase() {
      Rbayz::parList.push_back(&par);
   };

   virtual ~modelBase() {  }
   
   // sample and update methods, now updating of hyper pars (varcomps) is separated to
   // allow running with fixed hyper parameters.
   virtual void sample() = 0;
   virtual void sampleHpars() = 0;
   // virtual void updateGD() = 0; // updates using GD?
   virtual void restart() = 0;

   // prepForOutput is for model classes that need to make a transform of
   // parameters for output, the base class defines an 'empty' version.
   virtual void prepForOutput() { };

   parVector* par=0;

};

#endif /* modelBase_h */
