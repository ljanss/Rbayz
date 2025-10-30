//  R/bayz
//  dataFactor.h - classes to store factor data with possible interactions, from model terms like A:B:C.
//  There are two versions:
//   - dataFactor: derives from simpleFactor and recodes interactions so it is again presented as a 
//     regular factor using the simpleFactor setup but with new level coding and labels that have the
//     combinations of all levels in the different factors. 
//   - dataFactorNC (not collapsed version): does NOT derive from simplefFacor and does not recode
//     interactions, instead it keeps holding multiple factors in a vector<simpleFactor *>.
//     The rn_cor models without mergeKernels can work with this vector of factors.
//
// storing one or multiple interacting factors with coding of the interaction levels
//     and making merged labels like "A1.B1.C1". 
//     Derives from simpleFactor and can also store the individual factors in a vector<simpleFactor*>
//     factorList.
//     This class is a bit more complex because it can handle several cases in the background and make
//     it useable for different model types:
//      - for objects like modelFixf, the factorList is temporarily used, then interactions are recoced
//       in the main object member variables data, nelem, labels (member variables deriving from simpleFactor)
//      - model objects like model_rn_cor, can use the factorList, without recoded interactions,
//        then the main member variables like data, labels remain uninitialised!
//     The setup is therefore somewhat dangerous if not used consistently. Model object constructors can
//     toggle the collapsing with the collapseInteractions argument, and the object then sets 'collapsed'
//     true. In that case the main member variables are initialised and should be used; if collapsed==false,
//     the factorList should be used. Rbayz does not have a high level of protection around all the member
//     variables, so if this is not used consistently you can make the code crash.
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
