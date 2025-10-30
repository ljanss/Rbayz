//  R/bayz
//  dataFactor.h - classes to store factor data with possible interactions, from model terms like A:B:C.
//  There are two versions:
//   - dataFactor: derives from simpleFactor and recodes interactions so it is again presented as a 
//     regular factor using the simpleFactor setup but with new level coding and labels that have the
//     combinations of all levels in the different factors. 
//   - dataFactorNC (Not Collapsed / recoded version): does NOT derive from simplefFacor and does not
//     recode interactions, instead it keeps holding multiple factors in a vector<simpleFactor *>.
//     The rn_cor models without mergeKernels can work with this vector of factors.
// If a model term has a single variable only (no interaction), this is handle with the standard
// dataFactor.
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

class dataFactorNC {

   dataFactorNC(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames, 
               std::vector<varianceSpec> varlist);
   ~dataFactorNC();
   int Nvar;
   std::vector<simpleFactor *> factorList;

};

#endif /* dataFactor_h */
