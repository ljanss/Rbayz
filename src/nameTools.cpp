//
//  BayzR --- nameTools.cpp
//
//  Created by Luc Janss on 02/07/2020.
//

#include "nameTools.h"
#include "rbayzExceptions.h"
#include <algorithm>
#include <map>

//#include "dataMatrix.h"
//#include "dataFactor.h"

// Transform R CharacterVector to the c++ equivalent vector<string>
// Note: the C++ vector should be empty, the algorithm uses push_back.

void CharVec2cpp(std::vector<std::string> & CppStrings, Rcpp::CharacterVector RStrings) {
   std::string s;
   for(size_t i=0; i< RStrings.size(); i++) {
      s = Rcpp::as<std::string>(RStrings[i]);
      CppStrings.push_back(s);
   }
}

// Comparison function to search for a string in a vector of pair<string, number>
inline bool compString2Pair(const std::pair<std::string, size_t> & p, const std::string & s)
{
    return p.first < s;
}

// Get rownames from an R NumericMatrix, and return in a C++ vector of strings.
// Throw errors when dimnames is missing or dimnames[0] (the rownames) is NULL.

void getMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat) {
   Rcpp::CharacterVector matNames;
   if (mat.hasAttribute("dimnames")) {
      Rcpp::List dimnames = Rcpp::as<Rcpp::List>(mat.attr("dimnames"));
      if(dimnames[0]!=R_NilValue) {
         matNames = dimnames[0];
      }
      else if (dimnames[0]!=R_NilValue) {
            matNames = dimnames[1];
      }
      else {
         throw(generalRbayzError("There are no rownames or colnames on matrix input"));
      }
   }
   else {
      throw(generalRbayzError("There are no col and rownames on matrix input"));
   }
   if(mat.nrow() != matNames.size())
      throw(generalRbayzError("Dimnames corrupted: length does not match matrix row size"));
   CharVec2cpp(names, matNames);
   
}

// Find 'name' in the column-names of a data frame.
// The column names are not sorted, so applying a simple sequential search, but it will
// get slow if the data frame is large (e.g. with large covariate data in it!).
// The c++ find will do the same and does not return the index, so I just made a simple
// search myself that goes over the vector of names.
// Could be improved by making a sorted version of the column names and using binary_search,
// but ideally then also storing and re-using the sorted version for all look-ups.
// Then the sorted names must be prepared in the main functinon....?
int findDataColumn(Rcpp::DataFrame d, std::string name) {
   std::vector<std::string> colnames;
   CharVec2cpp(colnames, d.names());
   int col;
   for(col=0; col<colnames.size(); col++) {
      if(colnames[col]==name) break;
   }
   if(col==colnames.size())
      return -1;
   else
      return col;
}

void builObsIndex(std::vector<size_t> & obsIndex, dataFactor *F, dataMatrix *M) {
   int errors=0;
//   Rcpp::Rcout << "Going to resize obsIndex to " << F->data.size() << "\n";
   obsIndex.resize(F->data.size(),0);
   // Build a sorted list of the matrix labels paired with matrix entry-rows
   std::vector< std::pair<std::string, size_t> > matLabelsSorted;
   for(size_t i=0; i< M->labels.size(); i++)
      matLabelsSorted.push_back(std::make_pair(M->labels[i], i));
   std::sort(matLabelsSorted.begin(),matLabelsSorted.end()); // shoud be enough
   // Make the index that links observations to matrix rows
   std::string s;
   std::vector< std::pair<std::string, size_t> >::iterator it;
   for(size_t i=0; i<F->data.size(); i++) {
      s = F->labels[F->data[i]];
      if ( (it = std::lower_bound(matLabelsSorted.begin(), matLabelsSorted.end(), s, compString2Pair)) != matLabelsSorted.end()) {
         obsIndex[i]=it->second;
      }
      else {
         // error
         Rcpp::Rcout << "ID " << s << "is not available in similarity/kernel matrix\n";
         errors++;
      }
   }
   if(errors)
      throw (generalRbayzError("Some IDs not found in similarity/kernel matrix"));
}