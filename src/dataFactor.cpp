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
   }
   else { // multiple factors to collapse: recode interaction levels
      // first build vector of combined labels matching the data
      size_t Ndata = factorList[0]->nelem;
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
      initWith(Ndata, 0);
      for(size_t i=0; i<Ndata; i++) {  // code the data
         p = new_unique_labels.find(new_data_labels[i]);
         data[i] = p->second;
      }
      // fill labels vector
      labels.reserve(lev);
      for(p=new_unique_labels.begin(); p != new_unique_labels.end(); p++) labels.push_back(p->first);
   }

   // remove the allocated factorList objects! The factorList going out of scope will not
   // be enough to clean up the objects pointed to.
   for(size_t i=0; i<factorList.size(); i++)
      delete factorList[i];

}

dataFactor::~dataFactor() {
}

/* ------------------- dataFactorNC class ------------------- 

note: the dataFactorNC class has chunks of code identical to dataFactor; it was first intergrated in one
class, but that got messy with object constructors and with dataFactorNC not needing the simpleFactor set-up,
so it is better as two classes. To avoid the identical code blocks, parts could be moved into a function
used by both constructors ...

*/


dataFactorNC::dataFactorNC(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames, 
               std::vector<varianceSpec> varlist)
{

   // Note: dataFactorNC should not be used for a single variable (no interaction) - the right triage should be done
   // in rbayz main to send that to another constructor. Then the code for 'non collapsed' interations does not have
   // to worry about handling the case of one factor with no interaction.
   if(variableObjects.size==1) {
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

   // Load and code the individual factors in the factorList (here a member variable).
   for(size_t i=0; i<variableObjects.size(); i++) {
      // If a factor has an associated kernel, the levels from the kernel are used in coding the factor - this
      // prepares to predict levels in the kernel that are not present in the data.
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

}


dataFactorNC::~dataFactorNC() {
   for(size_t i=0; i<factorList.size(); i++)
      delete factorList[i];
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
