//  simpleFactor
//

#include "simpleFactor.h"
#include "nameTools.h"
#include "rbayzExceptions.h"

// ----------------- simpleFactor class --------------------

/* [ToDo]?
   Missing levels in a factor now can go through and become a separate "NA" level.
   This may not always make sense or be desired.
   In the version with levelLabels from a kernel, "NA" will go through if "NA" is
   also in the levelLabels (is a rowname of the kernel), otherwise it will be an error
   of "unmatched levels" because "NA" will not be found in levelLabels.
*/

// comparison function for lower_bound on pair<string,int> by first element, used in the
// second contrusctor with levelLabels.
bool search_key_in_pair(const std::pair<std::string, int> & p, const std::string & key) {
    return p.first < key;
}

simpleFactor::simpleFactor(Rcpp::RObject col, std::string inp_name) : simpleIntVector()
{

   /* simpleFactor can take different types of input and convert it to a 'factor', but
      processing and handling NAs differs somewhat dependent on input type. I always
      keep NA as the last level.
      - if input is R factor, it is directly convertable to integer vector that can be
        copied in the 'data' using initWith, but R starts coding from 1 so it needs to
        subtract 1 to code from 0, and NAs will be large negative number and are repaired later.
      - if input is an R character vector, converting to strings would convert NAs to string "NA",
        but I avoid that by building a map with the input strings that are not NA. Then later
        NA is added as last level, and it's level-codes are inserted as the last 'lev' value.
      - if input is an R integer vector I build a map of <int,int> so the sorting is nicer (as
        strings it would give the ugly sorting with "10" before "2"); treatment of NAs is like for
        a character vector.
      - logical vector input can be directly converted to interger and copied in 'data' using
        initWith, NAs will become large negative numbers and needs NA treatment like the first case.
        The IntegerVector will have 0 for false, 1 for true, so it is immediate in C base-0 coding.
   */
   name = inp_name;
   if (Rf_isFactor(col))
   {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      initWith(Rtempvec); // init of 'data' vector
      Rcpp::CharacterVector templabels = col.attr("levels");
      CharVec2cpp(labels, templabels);
      for (size_t row = 0; row < unsigned(Rtempvec.size()); row++)
      {
         data[row] -= 1;
      }
      // If there was "NA" in the data, add it as extra (last) level
      if (Rcpp::sum(missing) > 0)
      {
         labels.push_back("NA");
         size_t last_level = labels.size() - 1;
         for (size_t row = 0; row < unsigned(Rtempvec.size()); row++)
         {
            if (missing[row])
               data[row] = last_level;
         }
      }
   }
   else if (Rcpp::is<Rcpp::IntegerVector>(col) && !Rf_isMatrix(col))
   {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      std::map<int, int> unique_levels;
      for (int i = 0; i < Rtempvec.size(); i++)
      { // build map with non-NA values
         if (!missing[i])
            unique_levels[Rtempvec[i]];
      }
      std::map<int, int>::iterator p;
      size_t lev = 0; // Code the merged levels in the map
      for (p = unique_levels.begin(); p != unique_levels.end(); p++)
         p->second = lev++;
      initWith(Rtempvec.size(), 0.0l);
      for (int i = 0; i < Rtempvec.size(); i++)
      { // code the data
         if (missing[i])
            data[i] = lev;
         else
         {
            p = unique_levels.find(Rtempvec[i]);
            data[i] = p->second;
         }
      }
      // fill labels vector
      labels.reserve(lev);
      for (p = unique_levels.begin(); p != unique_levels.end(); p++)
         labels.push_back(std::to_string(p->first));
      if (sum(missing) > 0)
         labels.push_back("NA");
   }
   else if (Rcpp::is<Rcpp::CharacterVector>(col) && !Rf_isMatrix(col))
   {
      Rcpp::CharacterVector Rtempvec = Rcpp::as<Rcpp::CharacterVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      std::vector<std::string> Ctempvec;
      CharVec2cpp(Ctempvec, Rtempvec);
      std::map<std::string, int> unique_levels;
      for (size_t i = 0; i < Ctempvec.size(); i++)
      { // build map with non-NA values
         if (!missing[i])
            unique_levels[Ctempvec[i]];
      }
      std::map<std::string, int>::iterator p;
      size_t lev = 0; // Code the merged levels in the map
      for (p = unique_levels.begin(); p != unique_levels.end(); p++)
         p->second = lev++;
      initWith(Ctempvec.size(), 0);
      for (size_t i = 0; i < Ctempvec.size(); i++)
      { // code the data
         if (missing[i])
            data[i] = lev;
         else
         {
            p = unique_levels.find(Ctempvec[i]);
            data[i] = p->second;
         }
      }
      // fill labels vector
      labels.reserve(lev);
      for (p = unique_levels.begin(); p != unique_levels.end(); p++)
         labels.push_back(p->first);
      if (sum(missing) > 0)
         labels.push_back("NA");
   }
   else if (Rcpp::is<Rcpp::LogicalVector>(col) && !Rf_isMatrix(col))
   {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      initWith(Rtempvec); // init of 'data' vector
      labels.push_back("FALSE");
      labels.push_back("TRUE");
      if (Rcpp::sum(missing) > 0)
      {
         labels.push_back("NA");
         for (size_t row = 0; row < unsigned(Rtempvec.size()); row++)
         {
            if (missing[row])
               data[row] = 2;
         }
      }
   }
   else
   {
      throw generalRbayzError("Variable/data column is not convertable to a factor: " + name);
   }
}

// Constructor version with supplied level labels - so far the levels / labels coming from a kernel.
// This uses a bit different strategy to first get all factor data as strings, then code according to levelLabels.
// The coding will be in the order of levelLabels.
simpleFactor::simpleFactor(Rcpp::RObject col, std::string name, std::vector<std::string> levelLabels, 
      std::string kernel_name) : simpleIntVector()
   {

   // get the factor from the data as vector<string> for all input types (except boolean?)
   std::vector<std::string> temp_fac_strings;
   if (Rf_isFactor(col))
   {
      Rcpp::IntegerVector temp_fac_Rlevs = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(temp_fac_Rlevs);
      Rcpp::CharacterVector templabels = col.attr("levels");
      temp_fac_strings.resize(temp_fac_Rlevs.size());
      for (size_t row = 0; row < unsigned(temp_fac_Rlevs.size()); row++) {
         if (missing[row])
            temp_fac_strings[row] = "NA";
         else
            temp_fac_strings[row] = Rcpp::as<std::string>(templabels[temp_fac_Rlevs[row] - 1]); // R factor levels from 1! 
      }
   }
   // I think IntegerVector and CharacterVector can be treated the same way here because the
   // Rcpp::as<Rcpp:CharacterVector> will convert integers to strings automatically.
   else if ( (Rcpp::is<Rcpp::IntegerVector>(col) || Rcpp::is<Rcpp::CharacterVector>(col) ) && !Rf_isMatrix(col)) {
      Rcpp::CharacterVector Rcpp_strings = Rcpp::as<Rcpp::CharacterVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rcpp_strings);
      temp_fac_strings.resize(Rcpp_strings.size());
      for (size_t row = 0; row < unsigned(Rcpp_strings.size()); row++) {
         if (missing[row])
            temp_fac_strings[row] = "NA";
         else
            temp_fac_strings[row] = Rcpp::as<std::string>(Rcpp_strings[row]);
      }
   }
   else if (Rcpp::is<Rcpp::LogicalVector>(col) && !Rf_isMatrix(col)) {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      temp_fac_strings.resize(Rtempvec.size());
      for (size_t row = 0; row < unsigned(Rtempvec.size()); row++) {
         if (missing[row])
            temp_fac_strings[row] = "NA";
         else if (Rtempvec[row] == 0)
            temp_fac_strings[row] = "FALSE";
         else
            temp_fac_strings[row] = "TRUE";
      }
   }
   else {
      throw generalRbayzError("Variable/data column is not convertable to a factor: " + name);
   }

   // Now code the factor data according to the supplied levelLabels. Here instead of a map
   // use a sorted vector and lower_bound to find levels. The map was useful when expecting
   // (maybe many) duplicates, but the levelLabels should be unique already. Then it is easy
   // to just make a sorted version for searching.
   // Note: the final coding remains in the levelLabels order, not the sorted order.
   // Note: to find back the entry in levelLabels, a vector of pairs is needed which holds the
   // label and the original index. 
   std::vector< std::pair<std::string,int> > sorted_levels_pairs(levelLabels.size());
   for (size_t i = 0; i < levelLabels.size(); i++) {
      sorted_levels_pairs[i] = std::make_pair(levelLabels[i], i);
   }
   std::sort(sorted_levels_pairs.begin(), sorted_levels_pairs.end());
   std::vector<std::string> unmatched_levels; // to store any unmatched levels
   for(size_t i=0; i< temp_fac_strings.size(); i++) {
      std::vector< std::pair<std::string,int> >::iterator it;
      it = std::lower_bound(sorted_levels_pairs.begin(), sorted_levels_pairs.end(),
                            temp_fac_strings[i], search_key_in_pair);
      if( it == sorted_levels_pairs.end() || it->first != temp_fac_strings[i] ) {
         unmatched_levels.push_back( temp_fac_strings[i] );
      }
      else {
         data[i] = it->second;  // code according to original levelLabels order
      }
   }


   // if any unmatched levels, throw error
   if( unmatched_levels.size() > 0 ) {
      Rbayz::Messages.push_back("There are levels in factor " + name +
         " that cannot be matched to rownames of the kernel " + kernel_name + ":");
      size_t nshow = std::min( unmatched_levels.size(), 10 );
      std::string s;
      for(size_t i=0; i< nshow; i++) {
         s += unmatched_levels[i] + " ";
      }
      if ( unmatched_levels.size() > nshow ) {
         s += " [+ " + std::to_string(unmatched_levels.size() - nshow) + " more]";
      }
      Rbayz::Messages.push_back(s);
      Rbayz::needStop = true;
      throw generalRbayzError("Error matching kernel to factor levels - see messages output");
   }


}

// Convert stored factor data back to a vector of strings - same as R as.character(factor).
std::vector<std::string> simpleFactor::back2vecstring()
{
   std::vector<std::string> result(nelem);
   for (size_t i = 0; i < nelem; i++)
   {
      result[i] = labels[data[i]];
   }
   return result;
}