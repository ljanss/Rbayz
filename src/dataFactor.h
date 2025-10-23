//  R/bayz
//  dataFactor.h - storing one or multiple interacting factors with coding of the interaction levels
//     and making merged labels like "A1.B1.C1". 
//     Makes use of simpleFactor to code and store the single factors, the single factors are also kept
//     because some models make use of that (but for others it could be dropped ...).
//     If there is only one factor, dataFactor is basically a wrapper around one simpleFactor.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef dataFactor_h
#define dataFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "simpleFactor.h"
#include "optionsInfo.h"

class dataFactor {
public:
   /* Constructors with and without a variance-list, and for one Robject or vector<Robjects>;
      As far as I can see there is no need for the combination of one Robject and one variance-object,
      so that combination is not made ... */
   dataFactor(Rcpp::RObject variableObject, std::string variableName);
   dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames);
   dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames, 
               std::vector<varianceSpec> varlist);
   void run_constructor(std::vector<Rcpp::RObject> variableObjects, 
         std::vector<std::string> variableNames, std::vector<varianceSpec> varlist);
   ~dataFactor();
   simpleIntVector levcode;  // coded level info
   std::vector<std::string> labels;
   int Nvar;  // The number of variables (interactions) in this factor
   std::vector<simpleFactor *> factorList;
   simpleFactor* onefactor=0;
   size_t nelem;
};

#endif /* dataFactor_h */
