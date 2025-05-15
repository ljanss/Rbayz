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
#include "Rbayz.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

// a struct with option specifications
struct optionSpec
{
   std::string optionText; // the original option as text, used for error reporting
   bool isgiven=false;     // whether option is given; this is used to return a 'not given' option
   bool haserror=false;
   Rcpp::Robject varObject;
   int format;             // set to 0 as default for not yet defined
   std::string keyw="";
   std::string key2="";
   std::string valstring="";
   bool valbool=false;
   std::vector<double> valnumb;
   optionSpec() {          // default constructor is needed as soon as specifying additional constructors ... 
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

// a struct with variance specification
struct varianceSpec
{
   std::string optionText; // the original complete vartruct as text, used for error reporting
   bool haserror;
   bool iskernel;
   std::string keyw;
   std::vector<optionSpec> varOptions;
   varianceSpec() {
      optionText="";
      haserror=false;
      iskernel=false;
      keyw="";
   }
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
      {"MIXT","counts",true}
   };
   std::map<std::string, int> option2format
   {
      std::make_pair("trace",2),
      std::make_pair("save",2),
      std::make_pair("V",1),
      std::make_pair("prior",4),
      std::make_pair("vars",3),
      std::make_pair("counts",3)
   };
public:
   optionsInfo(std::string fname, std::string optstring);
   ~optionsInfo() {

   }
   optionSpec& operator[](std::string&); 
};


#endif /* optionsInfo_h */
