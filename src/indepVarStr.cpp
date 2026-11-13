#include "indepVarStr.h"

// ---- indepVarStr class ----

indepVarStr::indepVarStr(parsedModelTerm & modeldescr, parVector* cpar) 
                 : modelVar(modeldescr), weights() {
       coefpar = cpar;
       weights.initWith(coefpar->nelem, 1.0l);
}

// A generic variance sample / update method to estimate a scale by regressing
// fitted values on residuals. The estimate is still stored as variance in par->val[0].
// The lhs and rhs to be passed can be made by getFitScaleStats(lhs, rhs), which is a method
// for all modelCoeff classes.
// Note: at the start of an MCMC it can happen that lhs is zero (lhs is sum(fit^2)). One
// case is an MH sampler in mxi model where nothing has been accepted to change from the initial
// zero starting value, which can happen when initial variance (started at 1) is too large.
// I'm expecting this problem is mostly solved by now using a good initial variance in the
// concerned models (currently Rreg and Mixt), but if it still occurs, the code here make
// 0.8 factor shrinkages for a max of 20 times (about 0.01 of the original value) to try to adjust.
// [ToDo] There is also a possibilitiy the scale could get negative! This is now avoided as
// the stored variance is always scale^2 and positive, and the returned scale is abs(scale), but
// this may not be an appropriate MCMC update. A correct approach is probably to sample from
// truncated Normal. This is like saying it has a likelihood contribution that is a full Normal,
// but a prior that restricts it to be positive.
double indepVarStr::sampleScale(double lhs, double rhs) {
    double curr_scale = sqrt(par->val[0]);
    double scale;
    static int lhs_zero_adj=0;
    if(lhs == 0.0l) {
        if (lhs_zero_adj < 20) {
            scale = 0.8*curr_scale;
            lhs_zero_adj++;
        }
        else {
            throw(generalRbayzError("Problems getting the scale right in sampleScale"));
        }
        par->val[0] = scale*scale;
        return scale;
    }
    else {
       double sample_mean = curr_scale*rhs/lhs;
       double sample_sd = curr_scale/sqrt(lhs);
       scale = R::rnorm( sample_mean, sample_sd);
    }
    par->val[0] = scale*scale;
    double invvar = 1.0l/par->val[0];
    for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
    return abs(scale);
}

// ---- idenVarStr class ----

idenVarStr::idenVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
   par = new parVector(modeldescr, 1.0l, "var");
   par->traced=1;
   par->varianceStruct="IDEN";
}

idenVarStr::~idenVarStr() {
    delete par;
}

void idenVarStr::restart() {
   double invvar = 1.0l/par->val[0];
   for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
}

void idenVarStr::sample() {
  double ssq=0.0;
  for(size_t k=0; k < coefpar->nelem; k++)
     ssq += coefpar->val[k]*coefpar->val[k];
  par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
  double invvar = 1.0l/par->val[0];
  for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
}

// ---- diagVarStr class ----

/* diagVarStr is the model b~N(0,Ds^2), and D is diagonal, or on scalar level
   b_i ~ d_i s^2. It is the opposite of the weighted model where b_i ~ s^2 / w_i.
   diagVarStr is also internally used for the rn_cor models running regression on
   eigenvectors, and D has the eigenvalues. 
*/

// constructor with D passed as simpleDblVector (used internally in e.g. kernel model based on evecs)
diagVarStr::diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar, simpleDblVector & Ddiag)
        : indepVarStr(modeldescr, coefpar) {
    if(coefpar->nelem != Ddiag.nelem) {
        throw(generalRbayzError("ERROR dimension of DIAG does not fit random effect size"));
    }
    par = new parVector(modeldescr, 1.0l, "var");
    par->traced=1;
    par->varianceStruct="DIAG";
    // [ToDo]? this is copying the diagonal, but surely in the case when used for
    // eigenvector regression, this Ddiag are the eigenvalues that remain in the kernelMatrix
    // object, and it can be enough to 'wrap' the diag vector here around the other.
    // Probably, it can be assumed this is always OK when using this constructor where the
    // third arg is a simpleDblVector (then it must remain present somewhere else until cleanup).
    diag.initWith(Ddiag);
}

diagVarStr::~diagVarStr() {
    delete par;
}

// "regular" constructor that gets variance info from the parsed model description and D from an RObject
diagVarStr::diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar)
{
    if(modeldescr.varianceStruct!="DIAG")  // varianceStuct DIAG assures there is exactly one element in the
                                           // allOptions.Vlist(), and it is DIAG
        throw(generalRbayzError("Wrong call to diagVarStr with variance structure "+modeldescr.varianceStruct));
    
    std::vector<varianceSpec> varlist = modeldescr.allOptions.Vlist();
    // I still check again there is exactly one varlist element and that is has a 'varname' option
    if(varlist.size()!=1 || !varlist[0]["varname"].isgiven) {
        throw(generalRbayzError("ERROR V=DIAG[] is missing variable name in "+modeldescr.shortModelTerm));
    }
    try {
        dataCovar tempDiag(varlist[0].varOptions[0].varObject, false, false, varlist[0]["varname"].valstring);
        if(coefpar->nelem != tempDiag.nelem) {
            throw(generalRbayzError("ERROR dimension of DIAG does not fit random effect size"));
        }
        par = new parVector(modeldescr, 1.0l, "var");
        par->traced=1;
        par->varianceStruct="DIAG";
        diag.initWith(static_cast<const simpleDblVector&>(tempDiag));
    }
    catch(std::exception &err) {
        Rbayz::Messages.push_back(std::string(err.what()));
        throw generalRbayzError("Error occured in processing DIAG[" + varlist[0].varOptions[0].valstring + 
            "] - maybe not a numeric vector?");
        Rbayz::needStop=true;
    }
}

void diagVarStr::restart() {
   double invvar = 1.0l/par->val[0];
   for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
}

void diagVarStr::sample() {
  double ssq=0.0;
  for(size_t k=0; k < coefpar->nelem; k++)
     ssq += coefpar->val[k]*coefpar->val[k]/diag.data[k];
  par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
  double invvar = 1.0l/par->val[0];
  for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
}

// The scaleVarStr is currently quite limited, it sets-up a variance object with one parameter,
// but its sample() is disabled.
// It can be an idea to expand this to have a fixed 'grid' of proportional variances, which can
// be useful in mixture models and also informative for output. The variance in the output can
// then be scale*(grid), however it would produce a zero variance when one of the grid values is zero
// (which gives a bit ugly post-processing with missing postVar and MCMC stats).
scaleVarStr::scaleVarStr(parsedModelTerm & modeldescr, parVector* coefpar, double initvar)
        : indepVarStr(modeldescr, coefpar)
{
    par = new parVector(modeldescr, initvar, "var");
    par->traced=1;
    par->varianceStruct="scale";
}

scaleVarStr::~scaleVarStr() {
    delete par;
}

void scaleVarStr::restart() { }

// The class hierarchy requires defining sample(), but variance in this class needs to be
// updated using sampleScale(). The coefficient class sampleHpars() should do this appropriately.
void scaleVarStr::sample() {
    throw generalRbayzError("Incorrect calling of scaleVarStr::sample()");
}

/* ---- grid-LASSO ----
   The weights vector is not used, but no easy way to avoid it being allocated in the parent class.
   Maybe it can be used in future for a 'reweighted' LASSO version, supplying prior weights.
   The model par[0] is variance, but estimated as a scaling factor.
   The constructor stores pointers to residual and fit-vectors in the coefficient class that uses this
   class as variance structure.
*/

gridLVarStr::gridLVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
    par = new parVector(modeldescr, 1.0l, "var");
    par->traced=1;
    par->varianceStruct="grLASS";
}

gridLVarStr::~gridLVarStr() {
     delete par;
}

void gridLVarStr::restart() {  }

// The class hierarchy requires defining sample(), but variance in this class needs to be
// updated using sampleScale(). The coefficient class sampleHpars() should do this appropriately.
void gridLVarStr::sample() {
    throw generalRbayzError("Incorrect calling of gridLVarStr::sample()");
}

/* ---- lassVarStr ---- Bayesian (Power?) LASSO
        b_i ~ Exp(- 'rate' |b_i|^'pow' ).
   'rate' is estimated from the data, 'pow' is tunable "power" parameter with default value 0.8.
   The 'rate' has a default improper uniform prior, or a user-supplied gamma prior.
   OBS: lasso implementation from bayz does not match the generic interface of idenVarStr
   with a vector of variances that is supplied to the 'lower' model to make mixed-model updates.
   ---> not yet finished, still in doubt if this is useful and if it would be the Park-Casella model
        with an extra layer of variances (which fits the common indepVarStr interface), or the
        old bayz implementation.
*/

// "regular" constructor that gets variance info from the parsed model description
lassVarStr::lassVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar)
{
    // need to check how to get the power parameter from the parsed modeldescr
    par = new parVector(modeldescr, 1.0l, "rate");
    par->traced=1;
    par->varianceStruct="LASS";
}

lassVarStr::~lassVarStr() {
    delete par;
}

void lassVarStr::restart() {
   double invvar = 1.0l/par->val[0];
   for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
}

// this is still copy from diagVarStr, but it will probably look most like idenVarStr ...
void lassVarStr::sample() {
  double ssq=0.0;
  for(size_t k=0; k < coefpar->nelem; k++)
     ssq += coefpar->val[k]*coefpar->val[k]/diag.data[k];
  par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
  double invvar = 1.0l/par->val[0];
  for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
}

// ---- logLinVarStr ---- (to do)

loglinVarStr::loglinVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
    // add constructor
}

loglinVarStr::~loglinVarStr() {
//    ....
}


void loglinVarStr::sample() {
    // add sample implementation
}

