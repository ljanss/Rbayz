//
//  dataFactor.cpp
//

#include "dataFactor.h"
#include "nameTools.h"
#include "rbayzExceptions.h"

dataFactor::dataFactor(Rcpp::RObject oneVarObject, std::string OneVarName) :
      levcode() {
   std::vector<Rcpp::RObject> temp_var_objects = {oneVarObject};
   std::vector<std::string> temp_var_names = {OneVarName};
   run_constructor(temp_var_objects, temp_var_names);
}

dataFactor::dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames) :
      levcode() {
   run_constructor(variableObjects, variableNames);
}

void dataFactor::run_constructor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames) {
   if(variableObjects.size()==1) {                                         // for only one factor use 'onefactor' to allocate a
      onefactor = new simpleFactor(variableObjects[0],variableNames[0]);   // simpleFactor object and make levcode a wrapper.
      levcode.data = onefactor->data;                                      //  ... but it's tricky for clean-up (see destructor)
      levcode.nelem = onefactor->nelem;
      labels = onefactor->labels;                                          // this make a copy? - so here some double memory use
      Nvar = 1;
      nelem = levcode.nelem;
   }
   else {         // multiple interacting factors
      for(size_t i=0; i<variableObjects.size(); i++)
         factorList.push_back(new simpleFactor(variableObjects[i],variableNames[i]));
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
