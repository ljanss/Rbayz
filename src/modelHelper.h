//
//  Rbayz --- modelHelper.h
//
//  "helper" class to store additional modelling vectors.
//  There is always an association with a certain model-term and the main parameter in that model-term.
//  Thus this class using a modeldescr (from which basis of the parameter name can be extracted), and the
//  parameter vector after which the helper needs to be set up (with equal size and labels).
//  The helper parameter vector then gets a name as supplied namePrefix + same name as the 'fiendPar' supplied.
//
//  Created by Luc Janss on 27/03/2025.
//

#ifndef modelHelper_h
#define modelHelper_h

#include <Rcpp.h>
#include <unistd.h>
#include "modelBase.h"

class modelHelper : public modelBase {
   
public:

   modelHelper(parsedModelTerm & modeldescr, double initVal, parVector& friendPar, std::string namePrefix)
                 : modelBase() {
      par = new parVector(modeldescr, initVal, friendPar, namePrefix);
   }

   ~modelHelper() {
      delete par;
   }

   void sample() { }
   void sampleHpars() { }
   void restart() { }

};

#endif /* modelHelper_h */
