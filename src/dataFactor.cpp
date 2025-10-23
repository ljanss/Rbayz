//
//  dataFactor.cpp
//

#include "dataFactor.h"
#include "nameTools.h"
#include "rbayzExceptions.h"
#include "optionsInfo.h"
#include <map>

dataFactor::dataFactor(Rcpp::RObject oneVarObject, std::string OneVarName) :
      levcode() {
   std::vector<Rcpp::RObject> temp_var_objects = {oneVarObject};
   std::vector<std::string> temp_var_names = {OneVarName};
   std::vector<varianceSpec> temp_varlist(1); // empty varianceSpec object
   run_constructor(temp_var_objects, temp_var_names, temp_varlist);
}

dataFactor::dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames) :
      levcode() {
   std::vector<varianceSpec> temp_varlist(variableObjects.size()); // empty varianceSpec objects
   run_constructor(variableObjects, variableNames, temp_varlist);
}

dataFactor::dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames,
         std::vector<varianceSpec> varlist) : levcode() {
   run_constructor(variableObjects, variableNames, varlist);
}


// Note: run_constructor always gets a variance-list, but it can be empty objects. This will work fine,
// because empty varianceSpec objects have iskernel=false and kernObject=R_NilValue.
// The varianceSpec object is used to use the levels from a kernel for coding the factor if iskernel==true.

/* working notes: 
   - I think can remove the use of onefactor, that simplifies code to always use the factorList, also for
     one factor;
   - the use of 'onefactor' was to avoid storing the same levcode and labels once in the factor in the
     factorlist, and once in the 'overall' levcode and labels (used to code interactions for >1 factors).
     Modeling objects will use the 'overall' codes. Without 'onefactor', avoiding storing the same info twice
     can still be done with the same (dangerous) copying of pointers.
   - add new constructor for simpleFactor to use labels from a kernel
   - the real coding work to do for using the kernel labels will then be in simpleFactor
This was some old code for 1 factor:
      levcode.data = onefactor->data;
      levcode.nelem = onefactor->nelem;
      labels = onefactor->labels;  // this make a copy? - so here some double memory use
      Nvar = 1;
      nelem = levcode.nelem;
*/
void dataFactor::run_constructor(std::vector<Rcpp::RObject> variableObjects, 
          std::vector<std::string> variableNames, std::vector<varianceSpec> varlist) {
   bool canUseVarlist = true;
   if(variableObjects.size() != variableNames.size()) {
      throw generalRbayzError("Something wrong in building factor: size of objects and names do not match");
   }
   if(varlist.size() != variableObjects.size()) {
      // this is also checked in the model_rn_cor constructors, where there is more context for reporting
      // the error. Therefore here just ignore it, but set needStop and don't use the varlist.
      canUseVarlist = false;
      Rbayz::needStop = true;
   }
   for(size_t i=0; i<variableObjects.size(); i++) {
      // check varlist to see if labels from kernel should be used in coding factor
      factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
   }
   size_t Ndata=factorList[0]->nelem;
   for(size_t i=1; i<factorList.size(); i++) {  // double check that the sizes of the factors are identical
      if( factorList[i]->nelem != Ndata) {
         std::string s="Interacting factors do not have the same length:";
         for (size_t j=0; j<factorList.size(); j++) {
            s += " " + variableNames[j] + "(" + std::to_string(factorList[j]->nelem) + ")";
         }
         throw generalRbayzError(s);
      }
   }
   if(variableObjects.size()==1) {
      // can copy levcode and labels from factorList[0] in 'main' levcode and labels
   }
   else { // multiple factors: main levcode and labels are for the interaction
      // first build vector of combined labels matching the data
      std::vector<std::string> new_data_labels(factorList[0]->back2vecstring());
      for(size_t i=1; i<factorList.size(); i++) {
         std::vector<std::string> next_strings(factorList[i]->back2vecstring());
         for(size_t j=0; j<Ndata; j++)
            new_data_labels[j] += "." + next_strings[j];
      }
      // build map to find and code unique labels
      std::map<std::string, int> new_unique_labels;
      for(size_t i=0; i<new_data_labels.size(); i++)
            new_unique_labels[new_data_labels[i]];
      std::map<std::string, int>::iterator p;
      size_t lev=0;        // Code the merged levels in the map
      for(p=new_unique_labels.begin(); p != new_unique_labels.end(); p++) p->second = lev++;
      levcode.initWith(Ndata, 0);
      for(size_t i=0; i<Ndata; i++) {  // code the data
         p = new_unique_labels.find(new_data_labels[i]);
         levcode[i] = p->second;
      }
      // fill labels vector
      labels.reserve(lev);
      for(p=new_unique_labels.begin(); p != new_unique_labels.end(); p++) labels.push_back(p->first);
      Nvar=factorList.size();
      nelem=Ndata;
   }
   Nvar = factorList.size();
}

dataFactor::~dataFactor() {
   if(Nvar==1) {                   // levcode vector was used as 'wrapper'
      delete onefactor;
      levcode.nelem=0;             // this avoids triggering clean-up in levcode
   }
   else {
      for(size_t i=0; i<factorList.size(); i++)
         delete factorList[i];
   }
}

/*
void dataFactor::setupFirstVariable(Rcpp::RObject col) {
   // Here is place where conversion of IntegerVector and CharacterVector to
   // factor could be added. 
   if (!Rf_isFactor(col)) {
      throw generalRbayzError("Variable is not a factor (unfortunately cannot get the name here)\n");
   }
   Rcpp::IntegerVector tempvec = Rcpp::as<Rcpp::IntegerVector>(col);
   Rcpp::LogicalVector missing = Rcpp::is_na(tempvec);
   data.initWith(tempvec);
   if (Rcpp::sum(missing) > 0) {
      labels.push_back("NA");
      for(size_t row=0; row < unsigned(tempvec.size()); row++) {
         if(missing[row]) data[row] = 0;
      }
   }
   else {
      for(size_t row=0; row < unsigned(tempvec.size()); row++) {
         data[row] -= 1;
      }
   }
   Rcpp::CharacterVector templabels = col.attr("levels");
   CharVec2cpp(labels, templabels);
   Nvar=1;
}
// addVariables adds another variable in a factor
void dataFactor::addVariable(Rcpp::DataFrame &d, size_t col) {
   Rcpp::RObject Rcol = Rcpp::as<Rcpp::RObject>(d[col]);
   addVariable(Rcol);
}
void dataFactor::addVariable(std::string varname) {
   Rcpp::Environment Renv;
   Rcpp::RObject Rcol = Renv[varname];
   addVariable(Rcol);
}
void dataFactor::addVariable(Rcpp::RObject Rcol) {
   if(Nvar==0)
      setupFirstVariable(Rcol);
   else { // add (interact) another factor with already stored factor(s)
      dataFactor tempFact(Rcol); 
      std::vector<std::string> oldlabels(labels);
      size_t nLevel1=labels.size();
      size_t nLevel2=tempFact.labels.size();
      labels.resize(nLevel1 * nLevel2);
      for(size_t i=0; i<nLevel1; i++) {  // generate the new labels
         for(size_t j=0; j<nLevel2; j++) {
            labels[i*nLevel2+j] = oldlabels[i] + "%" + tempFact.labels[j];
         }
      }
      // Replace existing data-level-coding with codes to match the new interaction
      for(size_t i=0; i<data.nelem; i++) {
         data[i] = data[i]*nLevel2 + tempFact.data[i];
      }
      Nvar++;
   }
}
*/
