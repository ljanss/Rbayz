//
//  dataFactor.cpp
//

#include "dataFactor.h"
#include "nameTools.h"
#include "rbayzExceptions.h"
#include "optionsInfo.h"
#include <map>

/* ------------------- dataFactor class ------------------- */

dataFactor::dataFactor(Rcpp::RObject oneVarObject, std::string OneVarName) :
      simpleFactor() {
   std::vector<Rcpp::RObject> temp_var_objects = {oneVarObject};
   std::vector<std::string> temp_var_names = {OneVarName};
   std::vector<varianceSpec> temp_varlist(1); // empty varianceSpec object
   run_constructor(temp_var_objects, temp_var_names, temp_varlist);
}

dataFactor::dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames) :
      simpleFactor() {
   std::vector<varianceSpec> temp_varlist(variableObjects.size()); // empty varianceSpec objects
   run_constructor(variableObjects, variableNames, temp_varlist);
}

dataFactor::dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames,
         std::vector<varianceSpec> varlist) : simpleFactor() {
   run_constructor(variableObjects, variableNames, varlist);
}


// Note: run_constructor always gets a variance-list, but it can be empty objects. This will work fine,
// because empty varianceSpec objects have iskernel=false and kernObject=R_NilValue.
// The varianceSpec object is used to use the levels from a kernel for coding the factor if iskernel==true.
// Note: dataFactor is derived from simpleFactor, so it has nelem, data, labels, name members, and uses
// the 'empty' simpleFactor constructor, so nothing is filled for these member variables yet.
// The run_constructor makes a vector<simpleFactor *> and that is based on the simpleFactor constructors
// that can recode many types of R objects into factors and that fill the member variables.

void dataFactor::run_constructor(std::vector<Rcpp::RObject> variableObjects, 
          std::vector<std::string> variableNames, std::vector<varianceSpec> varlist)
{

   bool canUseVarlist = true;
   if(variableObjects.size() != variableNames.size()) {
      throw generalRbayzError("Something wrong in building factor: number of objects and names do not match");
   }
   if(varlist.size() != variableObjects.size()) {
      // this is also checked in the model_rn_cor constructors, where there is more context for reporting
      // the error. Therefore here just ignore it, but set needStop and don't use the varlist.
      canUseVarlist = false;
      Rbayz::needStop = true;
   }

   std::vector<simpleFactor *> factorList; // here factorList is local

   // Load and code the individual factors in the factorList.
   for(size_t i=0; i<variableObjects.size(); i++) {
      // If a factor has an associated kernel, the levels from the kernel are used in coding the factor
      if( canUseVarlist && varlist[i].iskernel ) {
         Rcpp::NumericMatrix temp_kernel = Rcpp::as<Rcpp::NumericMatrix>(varlist[i].kernObject);
         std::vector<std::string> temp_rownames = getMatrixNames(temp_kernel, 1);
         if(temp_rownames.size()>0) {
            factorList.push_back(new simpleFactor(variableObjects[i], variableNames[i], temp_rownames, varlist[i].keyw));
         }
         else { // ignore rownames? It will create an error later, but need to check if it will not cause problems in coding here.
            factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
            Rbayz::Messages.push_back("Warning: cannot retrieve row names for kernel [" + varlist[i].keyw + "]");
            Rbayz::needStop = true;
         }
      }
      else // simpler type of factor without kernel
         factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
   }
   Nvar = factorList.size();

   // double check that the row-sizes of the factors are identical
   size_t Ndata=factorList[0]->nelem;
   for(size_t i=1; i<factorList.size(); i++) {  
      if( factorList[i]->nelem != Ndata) {
         std::string s="Interacting factors do not have the same length:";
         for (size_t j=0; j<factorList.size(); j++) {
            s += " " + variableNames[j] + "(" + std::to_string(factorList[j]->nelem) + ")";
         }
         throw generalRbayzError(s);
      }
   }

   // recode and copy factor data in the member variables.
   // Note: the initWith versions from simpleIntVector are a bit limited, can only use
   // initWith(size, scalar-value), then copy contents in data, improve? [ToDo]?.
   if(variableObjects.size()==1) {
      // can copy levcode and labels from factorList[0] in 'main' data and labels.
      initWith(factorList[0]->nelem, 0);
      for(size_t i=0; i<nelem; i++)
         data[i] = factorList[0]->data[i];
      labels = factorList[0]->labels;  // instead of copying, it could also 'move' contents ...?
      name = factorList[0]->name;
   }
   else { // multiple factors to collapse: recode interaction levels by pasting labels together.
      std::vector<std::string> pasted_data_labels;
      paste_data_labels(pasted_data_labels, factorList);
      std::vector<std::string> unique_labels = pasted_data_labels;
      std::sort(unique_labels.begin(), unique_labels.end());
      std::vector< std::string>::iterator last = std::unique (unique_labels.begin(), unique_labels.end());
      unique_labels.erase(last, unique_labels.end());
      initWith(Ndata, 0);
      std::vector<std::string>::iterator p;
      // code the data by searching the pasted_data_labels in the unique_labels
      for(size_t i=0; i<Ndata; i++) {
         p = std::lower_bound(unique_labels.begin(), unique_labels.end(), pasted_data_labels[i]);
         data[i] = p - unique_labels.begin();
      }
      // fill labels vector
      labels.reserve(unique_labels.size());
      for(p=unique_labels.begin(); p != unique_labels.end(); p++) labels.push_back(*p);
      // The name of the interaction coming from the model-term is lost here, reconstruct it ...
      name = factorList[0]->name;
      for(size_t i=1; i< factorList.size(); i++) name += ":" + factorList[i]->name;
   }

   // remove the allocated factorList objects! The factorList going out of scope will not
   // be enough to clean up the objects pointed to.
   for(size_t i=0; i<factorList.size(); i++)
      delete factorList[i];

}

dataFactor::~dataFactor() {
}

/* ------------------- dataFactorNC class ------------------- */

// dataFactorNC also derives from simpleFactor, and has similar working, but keeps the factorList of simpleFactor objects.
// So there is some code duplication here, but it was too messy to integrate both in one class.

dataFactorNC::dataFactorNC(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames, 
               std::vector<varianceSpec> varlist) : simpleFactor(), firstOccurence()
{

   // dataFactorNC is not for a single variable (no interaction)!
   if(variableObjects.size()==1) {
      throw generalRbayzError("Error wrong calling of dataFactorNC, pls report to developers");
   }

   bool canUseVarlist = true;
   if(variableObjects.size() != variableNames.size()) {
      throw generalRbayzError("Something wrong in building factor: number of objects and names do not match");
   }
   if(varlist.size() != variableObjects.size()) {
      // this is also checked in the model_rn_cor constructors, where there is more context for reporting
      // the error. Therefore here just ignore it, but set needStop and don't use the varlist.
      canUseVarlist = false;
      Rbayz::needStop = true;
   }

   // Load and code the individual factors in the factorList (here a member variable) - code repeated from dataFactor.
   for(size_t i=0; i<variableObjects.size(); i++) {
      if( canUseVarlist && varlist[i].iskernel ) {
         Rcpp::NumericMatrix temp_kernel = Rcpp::as<Rcpp::NumericMatrix>(varlist[i].kernObject);
         std::vector<std::string> temp_rownames = getMatrixNames(temp_kernel, 1);
         if(temp_rownames.size()>0) {
            factorList.push_back(new simpleFactor(variableObjects[i], variableNames[i], temp_rownames, varlist[i].keyw));
         }
         else { // ignore rownames? It will create an error later, but need to check if it will not cause problems in coding here.
            factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
            Rbayz::Messages.push_back("Warning: cannot retrieve row names for kernel [" + varlist[i].keyw + "]");
            Rbayz::needStop = true;
         }
      }
      else // simpler type of factor without kernel
         factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
   }
   Nvar = factorList.size();

   // double check that the row-sizes of the factors are identical
   size_t Ndata=factorList[0]->nelem;
   for(size_t i=1; i<factorList.size(); i++) {  
      if( factorList[i]->nelem != Ndata) {
         std::string s="Interacting factors do not have the same length:";
         for (size_t j=0; j<factorList.size(); j++) {
            s += " " + variableNames[j] + "(" + std::to_string(factorList[j]->nelem) + ")";
         }
         throw generalRbayzError(s);
      }
   }
   nelem = Ndata;

   // code the interaction levels - code repeated from dataFactor, but here don't need to handle the single factor case.
   // Additional: also fill firstOccurence vector while coding, it needs a helper vector to track if a level is already
   // seen or not to mark if a data row has the first occurence of a level, or else it is a duplicate.

   std::vector<std::string> pasted_data_labels;
   paste_data_labels(pasted_data_labels, factorList);
   std::vector<std::string> unique_labels = pasted_data_labels;
   std::sort(unique_labels.begin(), unique_labels.end());
   std::vector< std::string>::iterator last = std::unique (unique_labels.begin(), unique_labels.end());
   unique_labels.erase(last, unique_labels.end());
   initWith(Ndata, 0);
   std::vector<std::string>::iterator p;
   firstOccurence.initWith(Ndata, 0); // initialize all to 0 ('false')
   simpleIntVector seen_levels(unique_labels.size()); // to track seen levels, it is initialized to 0!
   // code the data by searching the pasted_data_labels in the unique_labels
   for(size_t i=0; i<Ndata; i++) {
      p = std::lower_bound(unique_labels.begin(), unique_labels.end(), pasted_data_labels[i]);
      data[i] = p - unique_labels.begin();
      // track first occurrence
      if (seen_levels[data[i]] == 0) {
         firstOccurence[data[i]] = 1; // 'true' for first occurence
      }
      seen_levels[data[i]] = 1;
   }
   // fill labels vector
   labels.reserve(unique_labels.size());
   for(p=unique_labels.begin(); p != unique_labels.end(); p++) labels.push_back(*p);
   // The name of the interaction coming from the model-term is lost here, reconstruct it ...
   name = factorList[0]->name;
   for(size_t i=1; i< factorList.size(); i++) name += ":" + factorList[i]->name;

}


dataFactorNC::~dataFactorNC() {
   for(size_t i=0; i<factorList.size(); i++)
      delete factorList[i];
}

// paste_data_labels: given a vector of simpleFactor pointers, create pasted labels for all data rows;
// The function is designed to avoid making a copy of the pasted labels vector, if it would be returned,
// now the calling code should make an empty vector and pass it by reference to be filled.
void paste_data_labels(std::vector<std::string> & pasted_labels, const std::vector<simpleFactor *> & factorList) {
   size_t Ndata = factorList[0]->nelem;
   pasted_labels = factorList[0]->back2vecstring();
   for(size_t i=1; i<factorList.size(); i++) {
      std::vector<std::string> next_strings(factorList[i]->back2vecstring());
      for(size_t j=0; j<Ndata; j++)
         pasted_labels[j] += "." + next_strings[j];
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
