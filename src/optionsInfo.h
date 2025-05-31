//
//  optionsInfo.h
//
// A class to store, check and retrieve options on a model-term (including options in its variance-structures).

#ifndef optionsInfo_h
#define optionsInfo_h

#include <string>
#include <vector>
#include <map>
#include <Rcpp.h>
#include <utility>

// a class to specify one option allowing to store a boolean, text, vector<double> or
// an Robject.
class optionSpec {
public:
   std::string optionText;  // the original option as text, used for error reporting
   bool isgiven=false;      // whether option is given; this is used to return a 'not given' option
   bool haserror=false;
   Rcpp::RObject varObject; // not variance but variable object!
   int format;              // set to 0 as default for not yet defined
   std::string keyw="";
   std::string key2="";
   std::string valstring="";
   bool valbool=false;
   std::vector<double> valnumb;
   optionSpec() {          // default constructor
      isgiven=true;
      haserror=false;
      varObject = R_NilValue;
      format=0;
      optionText="";
      keyw=""; key2="";
      valstring="";
      valbool=false;
   }
   optionSpec(bool b) {
      isgiven=b;
      haserror=false;
      varObject = R_NilValue;
      format=0;
      optionText="";
      keyw=""; key2="";
      valstring="";
      valbool=false;
   }
};

// a class to specify one variance-structure
class varianceSpec {
public:
   std::string optionText; // the original complete vartruct as text, used for error reporting
   std::string keyw;
   bool haserror;
   bool iskernel;          // extend to: is kernel, is reserved structure, is number?
   Rcpp::RObject kernObject;
   std::vector<optionSpec> varOptions;
   varianceSpec() {
      optionText="";
      haserror=false;
      iskernel=false;
      keyw="";
      kernObject = R_NilValue;
   }
   optionSpec operator[](std::string s);
};

/* A struct to store a combination of model-term, option, and 'required' flag.
   This is used to build the table modterm2option as a vector<modOptPair>.
   The operator== is defined so that two of these structs
   are equal when both modterm and option are equal.
*/
struct modOptPair {
   std::string modterm;
   std::string option;
   bool required;
   modOptPair(std::string m, std::string o, bool r) {
      modterm=m; option=o; required=r;
   }
   modOptPair(std::string m, std::string o) {
      modterm=m; option=o; required=false;
   }
   bool operator==(const modOptPair& other) const {
      return (modterm==other.modterm && option==other.option);
   }
};

class optionsInfo
{
private:
   std::vector<optionSpec> optionList;        // stores all options on a model-term with V= as string (if available)
   std::vector<varianceSpec> varstructList;   // stores decomposed variance structure (if V is given and not V=~)
   std::vector<modOptPair> modterm2option {
      {"mn","trace",false},
      {"fx","trace",false},
      {"rn","trace",false},
      {"rr","trace",false},
      {"mn","save",false},
      {"fx","save",false},
      {"rn","save",false},
      {"rr","save",false},
      {"rn","V",false},
      {"rr","V",false},
      {"rn","prior",false},
      {"rr","prior",false},
      {"MIXT","vars",true},
      {"MIXT","counts",true},
      {"KERN","dim",false},
      {"KERN","dimp",false},
      {"rn","alpha_est",false},
      {"rn","alpha_save",false}
   };
   std::map<std::string, int> option2format
   {
      std::make_pair("trace",4),
      std::make_pair("save",4),
      std::make_pair("V",1),
      std::make_pair("prior",6),
      std::make_pair("vars",5),
      std::make_pair("counts",5),
      std::make_pair("dim",3),
      std::make_pair("dimp",3),
      std::make_pair("alpha_est",4),
      std::make_pair("alpha_save",4)
   };
public:
   optionsInfo() { }
   void constr(std::string fname, std::string optstring);
   ~optionsInfo() { }
   bool haserror=false;
   optionSpec operator[](std::string);
   std::vector<varianceSpec>& Vlist() {
      return this->varstructList;
   }
   // There is no way (yet...) to directly retrieve values for an option. Issues in providing that are
   // 1) options can be text, bool or vector<double>, so it would need 3 retrieval functions?
   // 2) it was difficult to think how to return if option is not present. Text could return "", vector<double>
   // could return a zero-length vector, but missing bool value?
   // Now thinking it is still useful to have member variable or function:
   //  - istrue(): returns if bool option is set true, but can only determine if false together with isgiven
   //  - text(): return text (valstring) or "" when not given or not appropriate
   //  - ?? numeric values: there is instances where there is only one value to return, other cases will need the
   //       whole vector. Using a missing value could help for the one-value cases.
   // This would allow to use expressions like allOptions["V"].text(), or allOptions["save"].istrue(), etc.
};


#endif /* optionsInfo_h */
