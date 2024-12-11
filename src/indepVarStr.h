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
    indepVarStr(parsedModelTerm & modeldescr, parVector* cpar) : modelVar(modeldescr), weights() {
       coefpar = cpar;
       weights.initWith(coefpar->nelem, 1.0l);
   }
   simpleDblVector weights;
};

class idenVarStr : public indepVarStr {

public:

    idenVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
       par = new parVector(modeldescr, 1.0l, "var");
       par->traced=1;
       par->varianceStruct="IDEN";
    }

    ~idenVarStr() {
        delete par;
    }

    void restart() {
       double invvar = 1.0l/par->val[0];
       for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
       Rcpp::Rcout << "In idenVarStr weights set to " << invvar << "\n";
    }

    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefpar->nelem; k++)
         ssq += coefpar->val[k]*coefpar->val[k];
      par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
      double invvar = 1.0l/par->val[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
    }
};

/* diagVarStr is the model b~N(0,Ds^2), and D is diagonal, or on scalar level
   b_i ~ d_i s^2. It is the opposite of the weighted model where b_i ~ s^2 / w_i.
   diagVarStr is internally used for the rn_cor models running regression on eigenvectors,
   and D has the eigenvalues. 
*/
class diagVarStr : public indepVarStr {

public:

    // constructor with simpleDblVector (used internally in e.g. kernel model based on evecs)
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar, simpleDblVector & Ddiag)
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

    // "regular" constructor that gets variance info from the parsed model description
    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar)
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

    void restart() {
       double invvar = 1.0l/par->val[0];
       for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
    }

    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefpar->nelem; k++)
         ssq += coefpar->val[k]*coefpar->val[k]/diag.data[k];
      par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
      double invvar = 1.0l/par->val[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
    }

    simpleDblVector diag;

};

/* lassVarStr is the Bayesian Power LASSO, on scalar level it is the model
        b_i ~ Exp(- 'rate' |b_i|^'pow' ).
   'rate' is estimated from the data, 'pow' is tunable "power" parameter with default value 0.8.
   The 'rate' has a default improper uniform prior, or a user-supplied gamma prior.
   OBS: lasso implementation from bayz does not match the generic interface of idenVarStr
   with a vector of variances that is supplied to the 'lower' model to make mixed-model updates.
*/

class lassVarStr : public indepVarStr {

public:

    // "regular" constructor that gets variance info from the parsed model description
    lassVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar)
    {
        // need to check how to get the power parameter from the parsed modeldescr
        par = new parVector(modeldescr, 1.0l, "rate");
        par->traced=1;
        par->varianceStruct="LASS";
    }

    void restart() {
       double invvar = 1.0l/par->val[0];
       for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
    }

    // this is still copy from diagVarStr, but it will probably look most like idenVarStr ...
    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefpar->nelem; k++)
         ssq += coefpar->val[k]*coefpar->val[k]/diag.data[k];
      par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
      double invvar = 1.0l/par->val[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
    }

    simpleDblVector diag;

};

// mixtVarStr is now standard 2-class mixture with pi0, pi1, v0, v1
class mixtVarStr : public indepVarStr {
public:
    mixtVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
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

    void restart() {
        // ??
    }

    void sample() {
        // add sample implementation
    }

    int Ncat;
    std::vector<double> Vars;
    std::vector<int> Counts;
    simpleDblVector diag;

//    getSetOptions(modeldescr.options["V"]); // not yet defined?
//    modelBVS *mixmod;
};

class loglinVarStr : public indepVarStr {
public:
    loglinVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
        // add constructor
    }
    void sample() {
        // add sample implementation
    }
};

/* for weighted variance the SSQ statistic will be multiplied by weight
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]*weights[k];
*/

#endif /* indepVarStr_h */
