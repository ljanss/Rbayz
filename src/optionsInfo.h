//
//  optionsInfo.h
//
// A class to store, check and retrieve options on a model-term.



#ifndef optionsInfo_h
#define optionsInfo_h

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "rbayzExceptions.h"

// a struct with option specifications
struct optionSpec
{
   bool isgiven;      // whether option is given by user?? maybe don't need this
   bool haserror;
   int format=0;      // 0: not defined; see onenote notebook
   std::string keyw;
   std::string key2;
   std::string valstring;
   bool valbool;
   std::vector<double> valnumb;
};

// a struct with variance speficications
struct varianceSpec
{
   bool haserror;
   bool iskernel;
   std::string keyw;
   std::vector<optionSpec> varOptions;
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
   bool operator==(const modOptPair& other) const {
      return (modterm==other.modterm && option==other.option);
   }
};

class optionsInfo
{
private:
   std::vector<optionSpec> optionList;
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
   optionsInfo(std::string, std::string);
   ~optionsInfo();
//   optionSpec& operator[](std::string&);
private:
};

optionsInfo::optionsInfo(std::string fname, std::string optstring)
{
   // split and store the options.
   // The options cannot simply be separated on commas, because there can be commas inside options
   // like: V=XX[a=1,b=2],W=YY[...]
   // Approach is therefore to scan character by character, check open & close brackets and compute
   // 'open_close_brack_balance'. A proper splitting comma is a comma where open_close_brack_balance==0. 
   // Options are then stored in a map<string, string>, split on first "=".
   size_t pos1,pos2,pos3,pos4,pos5;
   if (optstring!="") {
      pos1=-1;                          // start of first option, it will move up 1 at the start of the loop
      pos2=0;                           // will move to comma after first option, or end of string
      pos3=optstring.size()-1;          // position of last character in optString
      int open_close_brack_balance=0;
      std::string tmpstring;
      size_t pos4, pos5, tmpstring_len;
      do {
         pos1++;                        // for proper continuation: after processing an option,
                                        // pos1 will be left standing on the splitting comma.
         while( !(open_close_brack_balance==0 && optstring[pos2]==',') && pos2<pos3) {
            if(optstring[pos2]=='(' || optstring[pos2]=='[')
               open_close_brack_balance++;
            if(optstring[pos2]==')' || optstring[pos2]==']')
               open_close_brack_balance--;
            pos2++;
         }
         if(pos2==pos3)         // pos2 on the last character
            tmpstring=optstring.substr(pos1,(pos2-pos1+1));
         else                   // pos2 is after the last character (of piece to extract)
            tmpstring=optstring.substr(pos1,(pos2-pos1));
         /* tmpstring is an isolated option, it can have two formats (note: all spaces are removed):
              keyword=value
              keyword(value1,value2)
            In the initial parsing all is treated as text and stored in map 'options'. Multiple values
            from the second format are stored as a string of the comma-separated values.
            Note: keyword=value and keyword(value) are the same, they both end up in the map as options[keyword]=value.
         */
         pos4 = tmpstring.find('=');       // check option string for equal sign and open-parenth
         pos5 = tmpstring.find('(');
         tmpstring_len = tmpstring.size();
         if(pos4 != std::string::npos) {   // keyword=value
            options[tmpstring.substr(0,pos4)]=tmpstring.substr(pos4+1,tmpstring_len-pos4-1);
         }
         else if(pos5 != std::string::npos) {  // keyword(value1,value2)
            options[tmpstring.substr(0,pos5)]=tmpstring.substr(pos5+1,tmpstring_len-pos5-2);
         }
         else {
            throw generalRbayzError("Error: option [" + tmpstring + "] is not keyword=value or keyword(value1,value2) in " +
               shortModelTerm);
         }
         if(optstring[pos2]==',') {    // the while will continue for a next option
            pos1=pos2;
            pos2++;
         }
      }
      while (optstring[pos1]==',');
   }

}

optionsInfo::~optionsInfo()
{
}

/*
optionsInfoStruct& optionsInfo::operator[](std::string& s) {

}
*/

#endif /* optionsInfo_h */
