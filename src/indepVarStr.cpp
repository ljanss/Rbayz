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
void indepVarStr::sampleScale(double lhs, double rhs) {
    double curr_scale = sqrt(par->val[0]);
    double sample_mean = curr_scale*rhs/lhs;
    double sample_sd = curr_scale/sqrt(lhs);
    double scale = R::rnorm( sample_mean, sample_sd);
    par->val[0] = scale*scale;
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
   diagVarStr is internally used for the rn_cor models running regression on eigenvectors,
   and D has the eigenvalues. 
*/

// constructor with simpleDblVector (used internally in e.g. kernel model based on evecs)
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

// "regular" constructor that gets variance info from the parsed model description
diagVarStr::diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar)
{
    // for the moment only handling "simple" DIAG structure where varVariable[0] should have
    // name, and varObject[0] should have Robject with the diagonal info. 
    if(modeldescr.varianceStruct!="DIAG")
        throw(generalRbayzError("Wrong call to diagVarStr with variance structure "+modeldescr.varianceStruct));
    
    try {
        dataCovar tempDiag(modeldescr.varObject[0], false, false);
        if(coefpar->nelem != tempDiag.nelem) {
            throw(generalRbayzError("ERROR dimension of DIAG does not fit random effect size"));
        }
        par = new parVector(modeldescr, 1.0l, "var");
        par->traced=1;
        par->varianceStruct="DIAG";
        diag.initWith(tempDiag);
    }
    catch(std::exception &err) {
        Rbayz::Messages.push_back(std::string(err.what()));
        throw(generalRbayzError("Error occured in processing DIAG["+modeldescr.varVariable[0]+"]"));
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

// ---- mixtVarStr ---- (not yet finished)

// mixtVarStr is now standard 2-class mixture with pi0, pi1, v0, v1

mixtVarStr::mixtVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
    // some work to parse options in the MIXT[...] term; there should be 'vars' and 'counts' ...
    std::vector<std::string> split_options = splitString(modeldescr.varOption[0],",");
    std::string vars_string="",counts_string="";
    for(size_t i=0; i<split_options.size(); i++){  // check for correct syntax vars=c(...) and cut out part between (...)
        if(split_options[i].substr(0,7)=="vars=c(") vars_string=split_options[i].substr(7,(split_options[i].size()-8));
        if(split_options[i].substr(0,9)=="counts=c(") counts_string=split_options[i].substr(7,(split_options[i].size()-10));
    }
    if(vars_string=="" || counts_string=="") {
        throw(generalRbayzError("MIXT["+modeldescr.varOption[0]+"] is missing vars=c() or counts=c() or it is not well formatted/spelled"));
    }
    std::vector<std::string> vars_values_strings = splitString(vars_string,",");  // the single values split but still as strings
    std::vector<std::string> counts_values_strings = splitString(counts_string,",");
    if(vars_values_strings.size() != counts_values_strings.size()) {
        throw(generalRbayzError("In MIXT["+modeldescr.varOption[0]+"] number of elements in vars and counts are not equal"));
    }
    Ncat = vars_values_strings.size();
    Vars.resize(Ncat,0.0l);
    Counts.resize(Ncat,0);
    int total_counts=0;
    try {
        for(size_t i=0; i<Ncat; i++) {
            Vars[i]=std::stod(vars_values_strings[i]);
            Counts[i]=std::stoi(counts_values_strings[i]);
            total_counts += Counts[i];
        }
    }
    catch(std::exception &err) {
        Rbayz::Messages.push_back(std::string(err.what()));
        throw(generalRbayzError("Error in reading vars or counts values from MIXT["+modeldescr.varOption[0]+"]"));
    }
    std::vector<std::string> temp_labels = generateLabels("pi",Ncat);
    temp_labels.insert(temp_labels.begin(),"var");
    par = new parVector(modeldescr, 1.0l, temp_labels);
    for(size_t i=1; i<=Ncat; i++)          // The pi's are initialized from the prior counts
        par->val[i] = Counts[i]/total_counts;
    par->traced=1;
    par->varianceStruct="MIXT";
}

mixtVarStr::~mixtVarStr() {
    delete par;
}

void mixtVarStr::restart() {
    // ??
}

void mixtVarStr::sample() {
    // add sample implementation
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

