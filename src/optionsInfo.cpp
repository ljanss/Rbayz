#include "optionsInfo.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

// Check if options are allowed for a model-term or variance-structure; 
// - also resolve formats 2-3-4
// - also check if option set as boolean (and unknown in the option2format list) can be a variable name
// NOTE: check_options() is not a optionInfo class member, so need to pass modterm2option and option2format tables;
//       alternatively, as a member function it could have direct access to these tables.
int check_options(std::string fn, std::vector<optionSpec>& opts, std::vector<modOptPair>& modterm2option,
                    std::map<std::string, int>& option2format) {
   int errors=0;
   for(size_t opt=0; opt<opts.size(); opt++) {
      if(!opts[opt].haserror) {
         std::map<std::string, int>::iterator it1 = option2format.find(opts[opt].keyw);
         if(it1==option2format.end()) {           // Unfound option can be a variable name!
            if(opts[opt].format==1) {             // But this check will not find more complex A:B combinations...
               Rcpp::RObject varobj = getVariableObject(opts[opt].keyw);
               if(varobj==R_NilValue) {
                  Rbayz::Messages.push_back("Unrecognized option <"+opts[opt].keyw+"> (misspelled?) in <" + opts[opt].optionText + ">");
                  Rbayz::needStop = true;
                  opts[opt].haserror=true;
                  errors++;
               }
               else { // OK, relabel this option as a variable-name
                  opts[opt].varObject = varobj;
                  opts[opt].format = 10;
                  opts[opt].valbool = false;
                  opts[opt].valstring = opts[opt].keyw;  // the name gets stored in valstring and
                  opts[opt].keyw = "varname";            // the keyw becomes 'varname'
               }
            }
         }
         else {
            // the option is known, resolve 234 format ...
            if(opts[opt].format==234) opts[opt].format = it1->second;
            // and check if the combination of funcname (or variance struct name) and option is OK
            modOptPair tempMod2Opt(fn,opts[opt].keyw);
            std::vector<modOptPair>::iterator it2 = find(modterm2option.begin(),modterm2option.end(),tempMod2Opt);
            if(it2==modterm2option.end()) {
               Rbayz::Messages.push_back("Misplaced option <"+opts[opt].keyw+"> (not used here) in <" + opts[opt].optionText + ">");
               Rbayz::needStop = true;
               opts[opt].haserror=true;
               errors++;
            }
         }
      }
   }
   return errors;
}

// Parse values that are still stored in string to the 'bool' or 'numb' slots.
// If str2dbl conversion goes wrong, it will set needStop flag in the background, 
int parse_option_values(std::vector<optionSpec>& opts) {
   int errors=0;
   for(size_t opt=0; opt<opts.size(); opt++) {
      if(!opts[opt].haserror) {
         if (opts[opt].format==3) {
            opts[opt].valnumb.push_back(str2dbl(opts[opt].valstring, opts[opt].optionText));
         }
         else if (opts[opt].format==4) {
            if(opts[opt].valstring=="y" || opts[opt].valstring=="TRUE" || opts[opt].valstring=="T" || opts[opt].valstring=="1") {
               opts[opt].valbool = true;
            }
            else if (opts[opt].valstring=="n" || opts[opt].valstring=="FALSE" || opts[opt].valstring=="F" || opts[opt].valstring=="0") {
               opts[opt].valbool = false;
            }
            else {  // not a good boolean
               Rbayz::Messages.push_back("Boolean option <" + opts[opt].optionText + "> requires TRUE/FALSE y/n T/F or 1/0");
               Rbayz::needStop = true;
               opts[opt].haserror=true;
               errors++;
            }
         }
         else if ( opts[opt].format==5 || opts[opt].format==6) {
            std::vector<std::string> splits = splitString(opts[opt].valstring,",");
            for(size_t i=0; i<splits.size(); i++) {
               opts[opt].valnumb.push_back(str2dbl(splits[i], opts[opt].optionText));
            }
         }
      }
   }
   return errors;
}

void optionsInfo::constr(std::string fname, std::string optstring)
{
   // split options on comma (using splitStringNested that ignores nested commas)
   std::vector<std::string> optionsStrings = splitStringNested(optstring);

   // fill optionList vector with default/empty slots for each option
   optionList.resize(optionsStrings.size());

   // parse/store each option in the optionList. In the optionList variance options (V=) are stored as string.
   size_t errors=0;
   int varstruct_index=-1;
   for(size_t opt=0; opt<optionsStrings.size(); opt++) {
      optionList[opt].optionText = optionsStrings[opt];
      size_t equal,parenth,equalAfterParenth,star,optlen;
      equal = optionsStrings[opt].find('=');
      parenth = optionsStrings[opt].find_first_of("([");
      equalAfterParenth = (parenth==std::string::npos) ? std::string::npos : optionsStrings[opt].find('=',parenth);
      star = optionsStrings[opt].find("*");
      bool varstruct = (optionsStrings[opt].substr(0,2)=="V=");
      optlen = optionsStrings[opt].size();
      if(varstruct) {
         optionList[opt].keyw="V";
         optionList[opt].format=2;
         optionList[opt].valstring=optionsStrings[opt].substr(2,optlen-2);
         varstruct_index=opt;
      }
      else {
         if(star!=std::string::npos) {                                     // when not a varstruct, star is not allowed!
            Rbayz::Messages.push_back("Badly formatted option <"+optionsStrings[opt]+">: misplaced asterix(es)");
            Rbayz::needStop = true;
            optionList[opt].haserror=true;
            errors++;
         }
         else if (equal==std::string::npos && parenth==std::string::npos) {  // format 1, interpret as bool set true
            optionList[opt].format=1;
            optionList[opt].keyw=optionsStrings[opt];
            optionList[opt].valbool=true;
         }
         else if (equal!=std::string::npos && parenth==std::string::npos) {  // format 2-3-4 keyw=..., where '...' is
            optionList[opt].format=234;                                      // first stored as string
            optionList[opt].keyw=optionsStrings[opt].substr(0,equal);
            optionList[opt].valstring=optionsStrings[opt].substr(equal+1,optlen-equal-1);
         }
         else if (equal==std::string::npos && parenth!=std::string::npos) {  // format 5 with keyw(...), again '...'
            optionList[opt].format=5;                                        // first stored as string
            optionList[opt].keyw=optionsStrings[opt].substr(0,parenth);
            optionList[opt].valstring=optionsStrings[opt].substr(parenth+1,optlen-parenth-2);

         }
         else if (equal!=std::string::npos && parenth!=std::string::npos && equal < parenth &&
                    equalAfterParenth==std::string::npos) {                 // format 6, again values first go in string
            optionList[opt].format=6;
            optionList[opt].keyw=optionsStrings[opt].substr(0,equal);
            optionList[opt].key2=optionsStrings[opt].substr(equal+1,parenth-equal-1);
            optionList[opt].valstring=optionsStrings[opt].substr(parenth+1,optlen-parenth-2);
         }
         else {
            // other formats in the notebook are not handled now: formats 7,8 are not used (but could be extensions and
            // involve presence of equalAfterParenth), formats 9 and 11 are variances and are already handled in the first if,
            // and for the moment format 10 is not supported.
            Rbayz::Messages.push_back("Badly formatted option <"+optionsStrings[opt]+">: syntax not recognized");
            Rbayz::needStop = true;
            optionList[opt].haserror=true;
            errors++;
         }
      }
   } // end for handling first parsing model-term options

   // if a variance specification was present, and it does not start with ~, the 'variance structs' (elements of
   // the variance description separated by *) are separated and stored in slots in varstructList - and each
   // varstructList item has its own options. At this point variance-description stored as string in optionList.
   if(varstruct_index >= 0 && optionList[varstruct_index].valstring[0]!='~') {
      std::vector<std::string> varstructStrings = splitString(optionList[varstruct_index].valstring,"*");
      varstructList.resize(varstructStrings.size());
      size_t parenth, parenth2, equal2, optlen;
      for(size_t i=0; i<varstructStrings.size(); i++) {
         varstructList[i].optionText=varstructStrings[i];
         parenth = varstructStrings[i].find_first_of("([");
         if(parenth==std::string::npos) {                  // no parenth, varstruct without options
            varstructList[i].keyw=varstructStrings[i];     // (or a number!?)
         }
         else {                                            // varstruct with options within () or []
            varstructList[i].keyw=varstructStrings[i].substr(0,parenth);
            std::string optstring=varstructStrings[i].substr(parenth+1,(varstructStrings[i].size()-parenth-2));
            std::vector<std::string> optStrings = splitString(optstring,",");
            varstructList[i].varOptions.resize(optStrings.size());
            for(size_t j=0; j<optStrings.size(); j++) {                         // parse and store info from each
               equal2=optStrings[j].find('=');                                  // optStrings in the varOptions slots.
               parenth2=optStrings[j].find_first_of("([");
               optlen=optStrings[j].size();
               if (equal2==std::string::npos && parenth2==std::string::npos) {  // here same parsing as above, but
                  varstructList[i].varOptions[j].format=1;                      // options within varstructs can only
                  varstructList[i].varOptions[j].keyw=optStrings[j];            // be format 1, 234, or 5.
                  varstructList[i].varOptions[j].valbool=true;
               }
               else if (equal2!=std::string::npos && parenth2==std::string::npos) {
                  varstructList[i].varOptions[j].format=234;
                  varstructList[i].varOptions[j].keyw=optStrings[j].substr(0,equal2);
                  varstructList[i].varOptions[j].valstring=optStrings[j].substr(equal2+1,optlen-equal2-1);
               }
               else if (equal2==std::string::npos && parenth2!=std::string::npos) {
                  varstructList[i].varOptions[j].format=5;
                  varstructList[i].varOptions[j].keyw=optStrings[j].substr(0,parenth);
                  varstructList[i].varOptions[j].valstring=optStrings[j].substr(parenth+1,optlen-parenth-2);
               }
            }
         }
         // Determine if varstruct is a 'kernel' or reserved varstruct keyword; if kernel find RObject.
         // This could be extended to a third case where var-keyw is a number (to fix variances), then store it as 'fixed value'
         // in special slot in remove from varstructList ... ?
         if( (varstructList[i].keyw=="DIAG" || varstructList[i].keyw=="MIXT" || 
               varstructList[i].keyw=="LASS" || varstructList[i].keyw=="VCOV")) {
            varstructList[i].iskernel=false;
         }
         else {
            varstructList[i].iskernel=true;
            varstructList[i].kernObject = getVariableObject(varstructList[i].keyw);
            if(varstructList[i].kernObject == R_NilValue) {
               Rbayz::Messages.push_back("Variance structure <" + varstructList[i].keyw + "> is neither a known reserved structure nor a kernel matching an R object; misspelled?");
               Rbayz::needStop = true;
               varstructList[i].haserror=true;
               errors++;
               }
         }
      }
   }

   // check options in the main optionList and in each varstruct's varOptions list; this also resolves format 234
   // to be ready for the next step to parse values. The options are checked against the model-term funcName (fx, rn, etc.),
   // or for options in variance-structures against "KERN" if a kernel or the varstruct keyword otherwise (DIAG, MIXT, etc.)
   errors += check_options(fname, optionList, modterm2option, option2format);
   for(size_t i=0; i<varstructList.size(); i++) {
      if (varstructList[i].iskernel) errors += check_options("KERN", varstructList[i].varOptions, modterm2option, option2format);
      else errors += check_options(varstructList[i].keyw, varstructList[i].varOptions, modterm2option, option2format);
   }

   // Parse the values in all option-lists
   errors += parse_option_values(optionList);
   for(size_t i=0; i<varstructList.size(); i++) errors += parse_option_values(varstructList[i].varOptions);

   // Checking required options - now simply checking the few cases with required options and not using the required
   // flag in the modterm2option list. With only few requirements this looks easier, but clearly shows some repetitive coding ...
   for(size_t i=0; i<varstructList.size(); i++) {
      if(varstructList[i].keyw=="MIXT") {
         bool vars_present=false;
         bool counts_present=false;
         for(size_t j=0; j<varstructList[i].varOptions.size(); i++) {
            if (varstructList[i].varOptions[j].keyw=="vars") vars_present=true;
            if (varstructList[i].varOptions[j].keyw=="counts") counts_present=true;
         }
         if(!vars_present) {
            Rbayz::Messages.push_back("Variance structure <" + varstructList[i].optionText + "> is missing vars() specification");
            Rbayz::needStop = true;
            varstructList[i].haserror=true;
            errors++;
         }
         if(!counts_present) {
            Rbayz::Messages.push_back("Variance structure <" + varstructList[i].optionText + "> is missing counts() specification");
            Rbayz::needStop = true;
            varstructList[i].haserror=true;
            errors++;
         }
      }
      else if (varstructList[i].keyw=="DIAG") {
         bool varname_present=false;
         for(size_t j=0; j<varstructList[i].varOptions.size(); i++) {
            if (varstructList[i].varOptions[j].keyw=="varname") varname_present=true;
         }
         if(!varname_present) {
            Rbayz::Messages.push_back("Variance structure <" + varstructList[i].optionText + "> is missing a variable name (or it was not recognized as an R object)");
            Rbayz::needStop = true;
            varstructList[i].haserror=true;
            errors++;
         }
      }
   }

   if(errors) {
      this->haserror=true;
   }

} // end of construction function

// Retrieving options from the optionsInfo and varianceSpec objects;
// If requested option (string s) is not found, an optionSpec object with isgiven false is returned.
optionSpec optionsInfo::operator[](std::string s) {
   for(size_t opt=0; opt<optionList.size(); opt++) {
      if(optionList[opt].keyw==s) return optionList[opt];
   }
   return optionSpec(false);
}

optionSpec varianceSpec::operator[](std::string s) {
   for(size_t opt=0; opt<varOptions.size(); opt++) {
      if(varOptions[opt].keyw==s) return varOptions[opt];
   }
   return optionSpec(false);
}

