// modelHelper.h
// "helper" class to store additional modelling vectors.
// Using an additional model-object for this, the 'par' vector gets automatically
// added on the list of model parameter vectors.

#ifndef modelHelper_h
#define modelHelper_h

#include <Rcpp.h>
#include <unistd.h>
#include "modelBase.h"

class modelHelper : public modelBase {
   
public:

   // for now using constructor that needs labels and namePrefix, using standard default 0 to set-up the
   // par-vector.
   modelHelper(parsedModelTerm & modeldescr, Rcpp::CharacterVector& labels, std::string namePrefix)
                 : modelBase() {
      // maybe I need to make a new parVector constructor, for the lasso case there is no Rcpp:CharacterVector
      // of labels, but a vector<string>, and the combination with a namePrefix now does not exist.
      // It could also be convenient to use data from existing parVector to create a new one?
      par = new parVector(modeldescr, 0.0l, labels, namePrefix);
   }

   ~modelHelper() {
      delete par;
   }

};

#endif /* modelHelper_h */
