//
//  BayzR --- nameTools.h
//
//  Tools to retrieve, match names, make indexes, etc.
//  Created by Luc Janss on 02/07/2020.
//
#ifndef nameTools_h
#define nameTools_h

#include <vector>
#include <string>
#include <Rcpp.h>
#include "Rbayz.h"
#include "labeledMatrix.h"
#include "dataFactor.h"

void CharVec2cpp(std::vector<std::string> & labels, Rcpp::CharacterVector templabels);
std::vector<std::string> getMatrixNames(Rcpp::NumericMatrix & mat, int dim);
std::vector<std::string> generateLabels(std::string text, int n);
int findDataColumn(std::string name);

#endif /* nameTools_h */
