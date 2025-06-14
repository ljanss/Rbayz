//
//  BayzR --- nameTools.cpp
//
//  Created by Luc Janss on 02/07/2020.
//

#include "nameTools.h"
#include "rbayzExceptions.h"
#include <algorithm>
#include <map>

using Rsize_t = long int;

// Transform R CharacterVector to the c++ equivalent vector<string>
// Note: the C++ vector should be empty, the algorithm uses push_back;
// if not empty, elements will be added after existing elements!
void CharVec2cpp(std::vector<std::string> & CppStrings, Rcpp::CharacterVector RStrings) {
   CppStrings.reserve(CppStrings.size()+RStrings.size());
   for(Rsize_t i=0; i< RStrings.size(); i++) {
      CppStrings.push_back(Rcpp::as<std::string>(RStrings[i]));
   }
}

// getMatrixNames: attempts to retrieve row or col-names (dim=1 or 2) from an R matrix
// and return in an c++ vector<string>.
// Failure can be checked by the return vector to have size() 0.
std::vector<std::string> getMatrixNames(Rcpp::NumericMatrix & mat, int dim) {
   std::vector<std::string> names;
   if (mat.hasAttribute("dimnames")) {
      Rcpp::List dimnames = Rcpp::as<Rcpp::List>(mat.attr("dimnames"));
      if(dim==1 && dimnames[0] != R_NilValue) {
         Rcpp::CharacterVector matNames = Rcpp::as<Rcpp::CharacterVector>(dimnames[0]);
         CharVec2cpp(names, matNames);
      }
      else if (dim==2 && dimnames[1] != R_NilValue) {
         Rcpp::CharacterVector matNames = Rcpp::as<Rcpp::CharacterVector>(dimnames[1]);
         CharVec2cpp(names, matNames);
      }
   }
   // in all other cases (no dimnames, or dimnames empty), names is an empty vector
   return names;
}

// This is doing what R paste0(text,1:n) would do, i.e. generateLabels(kk,5) would generate kk1, kk2, ..., kk5
std::vector<std::string> generateLabels(std::string text, int n) {
   Rcpp::IntegerVector seq_ints = Rcpp::seq_len(n);      // 1..n as integers
   Rcpp::CharacterVector seq_strings                     // 1..n in Rcpp charvec
                   = Rcpp::as<Rcpp::CharacterVector>(seq_ints);
   std::vector<std::string> labels;
   for(int i=0; i<n; i++)                             // "text"+1..n
      labels.push_back(text+Rcpp::as<std::string>(seq_strings[i]));
   return labels;
}

// Find 'name' in the column-names of a data frame.
// The column names are not sorted, so applying a simple sequential search, but it will
// get slow if the data frame is large (e.g. with large covariate data in it!).
// The c++ find will do the same and does not return the index, so I just made a simple
// search myself that goes over the vector of names.
// Could be improved by making a sorted version of the column names and using binary_search,
// but ideally then also storing and re-using the sorted version for all look-ups.
// Then the sorted names must be prepared in the main function....?
int findDataColumn(std::string name) {
   std::vector<std::string> colnames;
   CharVec2cpp(colnames, Rbayz::mainData.names());
   size_t col;
   for(col=0; col<colnames.size(); col++) {
      if(colnames[col]==name) break;
   }
   if(col==colnames.size())
      return -1;
   else
      return col;
}
