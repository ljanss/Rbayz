//
//  parseFunctions
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef parseFunctions_h
#define parseFunctions_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>

void removeSpaces(std::string &s);
std::vector<std::string> splitString(std::string text, std::string splitchar);
std::vector<std::string> splitStringNested(std::string text);
int str2int(std::string s, std::string context);
double str2dbl(std::string s, std::string context);
std::string convertFormula(Rcpp::Formula f);
size_t findClosingBrack(std::string &s, size_t fromPos);
int getVariableType(Rcpp::RObject x);
Rcpp::RObject getVariableObject(std::string name);
std::vector<std::string> splitModelTerms(std::string mf);
//std::vector<std::string> getModelLHSTerms(std::string mf);
//std::vector<std::string> getModelRHSTerms(std::string mf);
std::vector<std::string> parseColNames(size_t col);
std::string getWrapName(std::string modelTerm);
std::string getFuncName(std::string modelTerm);
std::string getVarNames(std::string modelTerm);
std::string getVarDescr(std::string modelTerm);
std::string getPriorDescr(std::string modelTerm);
std::string getOptionText(std::string modelTerm, std::string text);
int checkOptions(std::map<std::string,std::string>, std::string);

#endif /* parseFunctions_h */
