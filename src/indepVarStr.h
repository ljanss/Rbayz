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
   double sampleScale(double lhs, double rhs);
   simpleDblVector weights;
   simpleDblVector& getWeights() { return weights; }
   // Now added querying methods to inform if the variance structure is IDEN
   // or DIAG. For these classes getScalarVar() returns the scalar variance,
   // and getDiagPtr() returns a pointer to the known diagonal part, or nullptr
   // for IDEN. For other classes getScalarVar() will return the first variance
   // in its parameter list (likely meaningless), and getDiagPtr() will return
   // nullptr.
   virtual bool isIDEN() const { return false; }
   virtual bool isDIAG() const { return false; }
   double getScalarVar() { return par->val[0]; }
   virtual const simpleDblVector* getDiagPtr() const { return nullptr; }
};

// The following two are currently in use for residual variance and fall in the group of
// 'simple' structures that are scalar * known diagonal. The main class that should accommodate
// more complex structures is the logLin structure (under development).

class idenVarStr : public indepVarStr {
public:
   bool isIDEN() const override { return true; }
public:
    idenVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~idenVarStr();
    void restart() override;
    void sample() override;
};

class diagVarStr : public indepVarStr {
public:
   bool isDIAG() const override { return true; }
   const simpleDblVector* getDiagPtr() const override { return &diag; }
public:
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar, simpleDblVector & Ddiag);
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~diagVarStr();
    void restart() override;
    void sample() override;
    simpleDblVector diag;
};

// scaleVarStr: a variance for random or shrinkage effects modelled as a scale.
// The variance structure has a single scale only, the real variance structure can be
// more complicated, but that is not relevant here.
// The scaleVarStr does not use the 'weight' vector used by other random effects (but
// it is difficult to avoid it being allocated in the parent class),
// and the standard sample() methods is disabled. Scale should be updated by calling
// the general updateScale() method with statistics that can be collected in coefficient
// model using getFitScaleStats() methods, but it may be possible to implement the
// update in the scaleVarStr.sample() to make it better organized and to allow more
// general treatment of priors.
// The scaleVarStr is currently implemted in modelMixt class, and the modelRanfi has
// test-code that has evaluated it on use for random effects.
// Note: scaleVarStr is not used / intended for the residual variance, and maybe can
// not be used for that because there is no 'fit' to regress on.

class scaleVarStr : public indepVarStr {
public:
    scaleVarStr(parsedModelTerm & modeldescr, parVector* coefpar, double initvar);
     ~scaleVarStr();
     void restart() override;
     void sample() override;
};

class gridLVarStr : public indepVarStr {
public:
    gridLVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
     ~gridLVarStr();
     void restart() override;
     void sample() override;
};

class lassVarStr : public indepVarStr {
public:
    lassVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~lassVarStr();
    void restart() override;
    void sample() override;
    simpleDblVector diag;
};

class loglinVarStr : public indepVarStr {
public:
    loglinVarStr(parsedModelTerm & modeldescr, parVector* coefpar);
    ~loglinVarStr();
    void sample() override;
};

#endif /* indepVarStr_h */
