//  parsedModelTerm.cpp

#include <Rcpp.h>
#include <vector>
#include <string>
#include "parsedModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

// parseModelTerm_step1: splits a model-term in 3 strings according to possible syntaxes
//  (1)   funcname(variableString,optionString)
//  (2)   funcname(variableString)
//  (3)   variableString
// Return is vector of 3 strings for funcname, variableString, optionString, some can
// be empty ("") if not available from syntaxes (2) and (3).
std::vector<std::string> parseModelTerm_step1(std::string mt) {

   std::vector<std::string> result(3);
   size_t pos1, pos2;
   bool has_funcname=false;

   // determine if there are parenthesis - so there should be a funcname
   if ( (pos1 = mt.find('(')) != std::string::npos)
      has_funcname=true;

   if (has_funcname && mt[mt.size()-1] != ')') {
      throw(generalRbayzError("No closing parenthesis in model-term: "+mt));
   }

   // get funcName and reset pos1 to where variable-list should start
   if ( has_funcname ) {
      result[0]=mt.substr(0, pos1);
      pos1++;
   }
   else {
      result[0]="";
      pos1=0;
   }

   // Determine if and where is ')' or ',' - this marks end of variableString with first 2
   // syntax patterns.
   // Checks: with a funcname, there must be ')' or ',' (but is already checked that there is ')')
   //         without a funcname, there cannot be ')' or ','
   pos2 = mt.find_first_of("),",pos1);
   if (!has_funcname && pos2 != std::string::npos) {     // wrong syntax
      if(mt[pos2]==')')
         throw(generalRbayzError("Unexpected closing parenthesis in response or model-term: "+mt));
      else
         throw(generalRbayzError("Unexpected comma in response or model-term: "+mt));
   }

   // 2. variableString
   size_t retrieve_length;
   if(pos2==std::string::npos)
      retrieve_length=pos2;            // 3rd syntax pattern: get substr to end
   else
      retrieve_length=pos2-pos1;       // 1st or 2nd syntax pattern: get substr to comma or )
   result[1]=mt.substr(pos1,retrieve_length);

   // 3. optionString: if pos2 is on a comma there are options upto closing parenthesis
   if(mt[pos2]==',') {
      retrieve_length = mt.size()-pos2-2;
      result[2]=mt.substr((pos2+1),retrieve_length);
   }
   else
      result[2]="";

   return result;

}

// parseModelTerm_step2: splitting / interpreting variables, options, etc.
// This one is defined as a member function to fill object member variables
void parsedModelTerm::parseModelTerm_step2(std::string fnName, std::string vrString, 
                     std::string optString, Rcpp::DataFrame &d) {

   funcName = fnName;
   variableString = vrString;

   // Make a shortened version of the model-term as text for messages
   if(vrString.length()<=12) {
     if(optString=="") shortModelTerm=funcName+"("+vrString+")";
     else shortModelTerm=funcName+"("+vrString+",...)";
   }
   else {
     shortModelTerm=funcName+"("+vrString.substr(0,12)+"...)";
   }

   // Analyse and split vrString
   size_t pos1 = vrString.find(':');
   size_t pos2 = vrString.find('/');
   size_t pos3 = vrString.find('|');
   if(pos2==std::string::npos && pos3==std::string::npos) {          // A, A:B, A:B:C
      if(pos1==std::string::npos) variablePattern="onevar";          // single variable, can be anything
      else variablePattern="intfactors";                             // A:B, A:B:C, must be interacting factors
   }
   else if (pos3!=std::string::npos && pos2==std::string::npos       // A|B, A|B:C
                  && (pos1==std::string::npos || pos2>pos3))         // but not A:B|C (not allowed)
      variablePattern="nestedreg";
   else if (pos2!=std::string::npos && pos1==std::string::npos       // A/B but no other patterns
                  && pos3==std::string::npos)                        // with | or : allowed with /
      variablePattern="rrcovars";
   else {
      std::string s="Cannot interpret/use variable specification \'"+vrString+
                    "\' in model-term: "+shortModelTerm;
      throw(generalRbayzError(s));
   }
   variableNames = splitString(vrString,":|/");

   // For every variable get an RObject pointing to it (whether it is from the data frame
   // or from R environment), and also store the type (factor, numeric, etc.) of the variable.
   for(size_t i=0; i<variableNames.size(); i++) {
      if (variableNames[i]=="1" || variableNames[i]=="0") {
         variableObjects.push_back(R_NilValue);
         variableTypes.push_back(0);
      }
      else {
         variableObjects.push_back(getVariableObject(d,variableNames[i]));
             // getVariableObject searches both the data frame 'd' and the R environment
         if(variableObjects.back() != R_NilValue)
            variableTypes.push_back(getVariableType(variableObjects.back()));
         else {
            throw generalRbayzError("Variable not found in data frame or R environment: "+variableNames[i]);
         }
      }
   }

   // split and store the options.
   // The options cannot simply be separated on commas, because there can be commas inside options
   // like: V=XX[a=1,b=2],W=YY[...]
   // Approach is therefore to scan character by character, check open & close brackets and compute
   // 'open_close_brack_balance'. A proper splitting comma is a comma where open_close_brack_balance==0. 
   // Options are then stored in a map<string, string>, split on first "=".
   if (optString!="") {
      pos1=-1;                          // start of first option, it will move up 1 at the start of the loop
      pos2=0;                           // will move to comma after first option, or end of string
      pos3=optString.size()-1;          // position of last character in optString
      int open_close_brack_balance=0;
      std::string tmpstring;
      size_t pos4;
      do {
         pos1++;                        // for proper continuation: after processing an option,
                                        // pos1 will be left standing on the splitting comma.
         while( !(open_close_brack_balance==0 && optString[pos2]==',') && pos2<pos3) {
            if(optString[pos2]=='(' || optString[pos2]=='[')
               open_close_brack_balance++;
            if(optString[pos2]==')' || optString[pos2]==']')
               open_close_brack_balance--;
            pos2++;
         }
         if(pos2==pos3)         // pos2 on the last character
            tmpstring=optString.substr(pos1,(pos2-pos1+1));
         else                   // pos2 is after the last character (of piece to extract)
            tmpstring=optString.substr(pos1,(pos2-pos1));
         // [ToDo] here thinking to also allow option syntax as "function-like" (looks like function),
         // which is xxxx(value1,value2). This is convenient to write a list of values.
         pos4 = tmpstring.find("=");  // locate equal sign in option string
         if(pos4 == std::string::npos) {
            throw generalRbayzError("Error: option [" + tmpstring + "] is not <keyword>=<value> in " + shortModelTerm);
         }
         options[tmpstring.substr(0,pos4)]=tmpstring.substr(pos4+1,std::string::npos);
         if(optString[pos2]==',') {    // the while will continue for a next option
            pos1=pos2;
            pos2++;
         }
      }
      while (optString[pos1]==',');
   }

   // Split and analyse the variance description. This writes in varianceStruct a string
   // that allows to select the right object class in main. Variances descriptions that are sequence of
   // variance-structures and kernels get split and annotated further in varType etc. 
   std::string variance_text = options["V"];  // also need to handle VE?
   if (variance_text=="") {
      varianceStruct="notgiven";
   }
   else {
      if (variance_text[0]=='~') {
         varianceStruct="llin";
      }
      else {    // all other cases should be structure keywords and kernels separated by stars
         std::vector<std::string> varianceElements = splitString(variance_text,"*");
         // Here there could be a way to allow fixing variances by detecting if the last
         // element is a numerical value, then store and remove that last element.
         // Split every variance element in a name and a parameter-part
         for(size_t i=0; i<varianceElements.size(); i++) {
            size_t bracket = varianceElements[i].find_first_of("([");
            size_t len_tot = varianceElements[i].length();
            std::string name, variable, options;
            if (bracket == std::string::npos) {   // simple variance-term like "Gmat"
               name = varianceElements[i];
               options = "";
               variable = "";
            }
            else {                                // variance-term with [...] like K1[dim=5] or DIAG[W]
               size_t closeBrack = findClosingBrack(varianceElements[i], bracket);
               if(closeBrack != (len_tot-1) ) {
                  throw generalRbayzError("Unbalanced parentheses in: "+varianceElements[i]);
               }
               name = varianceElements[i].substr(0,bracket);
               options = varianceElements[i].substr(bracket+1,(len_tot-bracket-2));
               if(name=="DIAG") {                 // separating first element in [...] as the variable for DIAG structures
                  size_t comma = options.find(',');     // [ToDo]to be extended to other structures where the first element
                  if(comma==std::string::npos) {        // is expected to be a variable name/object ...
                     variable = options;
                     options = "";
                  }
                  else {
                     variable = options.substr(0,comma);
                     options = options.substr(comma+1,std::string::npos);
                  }
               }
               else {                  // for not DIAG (and with extensions not other structures with a variable name/object),
                  variable="";         // variable is empty, options is all text between [... ] (possibly list of several comma-separated options) 
               }
            }
            // [ToDo] Here make varOption split in a elements (and further in map?). Now it is just the whole string ...
            varName.push_back(name);
            varVariable.push_back(variable);
            varOption.push_back(options);
         }  // end for-loop over varianceElements
         // if allowed to have DIAG in multiple var-elements, DIAG is not necessarily the last one,
         // and this check may need to move inside the above for-loop.
         if(varName.back()=="DIAG" && varVariable.back()=="") {
            throw generalRbayzError("DIAG specification is missing a variable in " + shortModelTerm);
         }
         // Fill the varianceObjects and varianceType vectors.
         // The Type is one of keywords IDEN, VCOV etc. OR "kernel";
         // The Objects are the kernel (matrix) object for type "kernel", the variable object for type "DIAG",
         // or R_NilValue in other cases. If an Robject is needed but cannot be found, that's an error.
         size_t nKernels=0, nVCOV=0;  // [ToDo] extend to also count the others ...
         for(size_t i=0; i<varianceElements.size(); i++) {
            std::string name=varName[i];
            if(name=="IDEN" || name=="MIXT" || name=="VCOV" || name=="WGHT" || name=="LLIN") {  
               varObject.push_back(R_NilValue);                            // I think LLIN cannot come in here ...
               varType.push_back(name);
               if(name=="VCOV") nVCOV++;
            }
            else if (name=="DIAG") {
               varObject.push_back(getVariableObject(d,varVariable[i]));
               if(varObject.back() == R_NilValue) {
                  throw generalRbayzError("Variable <" + varVariable[i]+">in DIAG[] not found in the R Environment in "
                                   + shortModelTerm);
               }
               varType.push_back(name);
            }
            else {   // expect a "kernel"
               varObject.push_back(getVariableObject(d,name));
               if(varObject.back() == R_NilValue) {
                  throw generalRbayzError("Variance/kernel <"+name+"> not found in the R Environment in "
                                   + shortModelTerm);
               }
               varType.push_back("kernel");
               nKernels++;
            }
         }
         // Determine the combination of variance-structures and set varianceStruct to indicate what
         // variance class to use for the modelling [ToDo] ... more work to be done here ...
         size_t nVarparts=varianceElements.size();
         if(nVarparts==1) {
            if (nKernels==1) varianceStruct="1kernel";
            if (nVCOV==1) varianceStruct="1VCOV";
            if (nKernels==0 && nVCOV==0) varianceStruct=varType[0];
         }
         else {  // multiple VTERMs
            if(nKernels==nVarparts)
               varianceStruct="kernels";
            else if (nKernels==(nVarparts-1) && nVCOV==1)
               varianceStruct="kernels-1vcov";
            else
               throw generalRbayzError("Cannot handle this variance pattern: "+variance_text);
         }
         // possible other cases to consider:
         // nkernels>0 && nVCOV>0 && (nKernels+nVCOV)==nVarparts:  kernels and multiple VCOV structures
         // cases with nKernels==0 and MIXT
         // cases with IDEN?
         // ...
      }  // end else (not llin var)
   } // end else (variance not empty)
} // end parsedModelTerm

// constructor for handling response term with separate variance description
parsedModelTerm::parsedModelTerm(std::string mt, std::string VEdescr, Rcpp::DataFrame &d)
{
   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);
   // for now not accepting functions on response, but it could be extended here to
   // allow e.g. log(Y), probit(Y), etc
   if(parse_step1[0]!="") throw(generalRbayzError("Unexpected function on response term: "+mt));
   // [workHere] the next one could just be a message
   if(parse_step1[2]!="") throw(generalRbayzError("Unexpected options retrieved from response term: "+mt));
   // here not inserted "rp" as funcName and residual variance description as ony possible option
   parseModelTerm_step2("rp", parse_step1[1], VEdescr, d);
}

// constructor for handling RHS model terms
parsedModelTerm::parsedModelTerm(std::string mt, Rcpp::DataFrame &d)
{
   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);
   parseModelTerm_step2(parse_step1[0], parse_step1[1], parse_step1[2], d);
}

std::ostream& operator<<(std::ostream& os, parsedModelTerm& p)
{
    os << p.shortModelTerm << ": funcName" << "[" << p.funcName << "] ";
    os << "variableString" << "[" << p.variableString << "] ";
    os << "variance" << "[" << p.options["V"] << "] ";
    os << "varianceStruct" << "[" << p.varianceStruct << "] ";
    os << "prior" << "[" << p.options["prior"] << "] ";
    os << "\n";
    return os;
}

/* variables that are in parsedModelTerm:
   std::string funcName="";
   std::string shortModelTerm="";
   std::string variableString="";
   std::string variablePattern="";
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   std::string varianceStruct="";
   std::string varianceLinMod="";
   std::vector<std::string> varOption;
   std::vector<std::string> varName;
   std::vector<Rcpp::RObject> varObject;
   std::vector<int> varianceKernelType;
   int hierarchType; // 0=no, 1=simplified form index/matrix, 2=genuine
   std::string hierarchModel="";
   std::string logging="";
*/