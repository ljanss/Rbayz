//  R/bayz
//  dataFactor.h - classes to store factor data with possible interactions, from model terms like A:B:C.
//  There are two classes, dataFactor and dataFactorNC which have subtle differences,
//  and are intended for different use by different model classes:
//   - dataFactor: recodes interactions and presents it as a regular single factor with new level coding and
//     labels; it does not keep the list of individual factors once the recoding is done. Used in Fixf, Ranfi,
//     Ranfc1 model classes.
//   - dataFactorNC (Not Collapsed version): this one keeps storing multiple factors in a vector<simpleFactor *>,
//     but also stores the recoded 'single factor' information, as well as a 'firstOccurence' vector. This
//     one is used by model classes working on the vector of factors, the single factor and firstOccurence info
//     is used in the backtransform. dataFactorNC is not intended for use for a single factor.
// Both versions can handle using levels from a kernel or multiple kernels where the factor is coded
// according to the kernel levels - this prepares to predict levels in the kernel that are not present
// in the data. This will also keep the levels ordered as they were given in the kernel.
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

class dataFactor : public simpleFactor {
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
   int Nvar;  // The number of variables (interactions) in this factor
};

class dataFactorNC : public simpleFactor {
public:
   dataFactorNC(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames, 
               std::vector<varianceSpec> varlist);
   ~dataFactorNC();
   int Nvar;
   std::vector<simpleFactor *> factorList;
   simpleIntVector firstOccurence; // similar as R duplicate() function, but opposite interpretation.
};

void paste_data_labels(std::vector<std::string> & pasted_labels, const std::vector<simpleFactor *> & factorList);

#endif /* dataFactor_h */
