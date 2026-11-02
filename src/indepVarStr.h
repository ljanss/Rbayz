// indepVarStr.h --- Classes of independent variance structures.
// These have a common "interface" of a vector of individual variances, which is
// defined in the base class and accessible through a pointer to base class.
// The actual form of this variance structure and implementation details for updating
// the variances can remain hidden for most (or all) uses of it.
// The log-linear variance model is also in this class, because it can interface in the
// same way with a vector of variances per random effect.

#ifndef indepVarStr_h
#define indepVarStr_h

#include <Rcpp.h>
#include "Rbayz.h"
#include "modelVar.h"
#include "modelBase.h"
#include "parsedModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"
#include "parVector.h"
#include "dataCovar.h"
#include "nameTools.h"
#include <unistd.h>

class indepVarStr : public modelVar {
public:
   indepVarStr(parsedModelTerm & modeldescr, parVector* cpar);
   virtual ~indepVarStr() { }
   void sampleScale(double lhs, double rhs);
   simpleDblVector weights;
};

class idenVarStr : public indepVarStr {
public:
    idenVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~idenVarStr();
    void restart();
    void sample();
};

class diagVarStr : public indepVarStr {
public:
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar, simpleDblVector & Ddiag);
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~diagVarStr();
    void restart();
    void sample();
    simpleDblVector diag;
};

// scaleVarStr: a variance for random or shrinkage effects modelled as a scale.
// The variance structure has a single scale only, the real variance structure can be
// more complicated, but that is not relevant here.
// The scaleVarStr does not use the 'weight' vector used by other random effects (but
// it is difficult to avoid it being allocated in the parent class),
// and the standard sample() methods is disabled. Scale should be updated by calling
// the general updateScale() method with statistics that can be collected in coefficient
// model using getFitScaleStats() methods. I made an implementation for the modelRregGRL
// class in modelRreg, but I will probably revise that one.

class scaleVarStr : public indepVarStr {
public:
    scaleVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
     ~scaleVarStr();
     void restart();
     void sample();
};

class gridLVarStr : public indepVarStr {
public:
    gridLVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
     ~gridLVarStr();
     void restart();
     void sample();
};

class lassVarStr : public indepVarStr {
public:
    lassVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~lassVarStr();
    void restart();
    void sample();
    simpleDblVector diag;
};

class mixtVarStr : public indepVarStr {
public:
    mixtVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~mixtVarStr();
    void restart();
    void sample();
    int Ncat;
    std::vector<double> Vars;
    std::vector<int> Counts;
    simpleDblVector diag;
//    getSetOptions(modeldescr.allOptions["V"]); // not yet defined?
//    modelBVS *mixmod;
};

class loglinVarStr : public indepVarStr {
public:
    loglinVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~loglinVarStr();
    void sample();
};

#endif /* indepVarStr_h */
