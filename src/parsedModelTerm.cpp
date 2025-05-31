//  parsedModelTerm.cpp

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

   std::vector<std::string> result(3,"");
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
   if(pos2 != std::string::npos && mt[pos2]==',') {
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
                     std::string optString) {

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
         variableObjects.push_back(getVariableObject(variableNames[i]));
             // getVariableObject searches both the data frame 'd' and the R environment
         if(variableObjects.back() != R_NilValue)
            variableTypes.push_back(getVariableType(variableObjects.back()));
         else {
            throw generalRbayzError("Variable not found in data frame or R environment: "+variableNames[i]);
         }
      }
   }

   allOptions.constr(funcName, optString);  // what if optString is empty, will all work OK?
   if(allOptions.haserror) {
      throw generalRbayzError("Errors in interpreting options in model-term " + shortModelTerm);
   }
   // Analyse the (combination of) variance structures. This writes in varianceStruct a string
   // that allows to select the right object class in main.
   optionSpec var_option = allOptions["V"];
   if( !var_option.isgiven ) {
      varianceStruct="notgiven";
   }
   else {
      std::string variance_text = var_option.valstring;
      if (variance_text[0]=='~') {
         varianceStruct="llin";
      }
      else {  // variance descriptions that are one or more variance-structures ...
         std::vector<varianceSpec> varianceList = allOptions.Vlist();
         size_t nKernels=0, nVCOV=0;                       // [ToDo] extend to also count the others ...
         for(size_t i=0; i<varianceList.size(); i++) {
            std::string name=varianceList[i].keyw;
            if(varianceList[i].iskernel) {
               nKernels++;
            }
            else {
               if (varianceList[i].keyw=="VCOV") {
                  nVCOV++;
               }
            }
         }
         size_t nVarparts=varianceList.size();
         if(nVarparts==1) {
            if (nKernels==1) varianceStruct="1kernel";
            if (nVCOV==1) varianceStruct="1VCOV";
            if (nKernels==0 && nVCOV==0) varianceStruct=varianceList[0].keyw;
         }
         else {  // multiple VTERMs - this now only accepts multiple kernels and is prepared to
            if(nKernels==nVarparts)                          // accept combinations with one VCOV 
               varianceStruct="kernels";
            else if (nKernels==(nVarparts-1) && nVCOV==1)
               varianceStruct="kernels-1vcov";
            else  // something mixed e.g. with reserved keyword structures
               varianceStruct="mixed";
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
parsedModelTerm::parsedModelTerm(std::string mt, std::string VEdescr)
{
   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);
   // for now not accepting functions on response, but it could be extended here to
   // allow e.g. log(Y), probit(Y), etc
   if(parse_step1[0]!="") throw generalRbayzError("Unexpected function on response term "+mt+" :"+parse_step1[0]);
   // [ToDo] the next one could just be a message
   if(parse_step1[2]!="") throw generalRbayzError("Unexpected options retrieved for response term "+mt+" :"+parse_step1[2]);
   // here not inserted "rp" as funcName and residual variance description as ony possible option
   parseModelTerm_step2("rp", parse_step1[1], VEdescr);
}

// constructor for handling RHS model terms
parsedModelTerm::parsedModelTerm(std::string mt)
{
   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);
   parseModelTerm_step2(parse_step1[0], parse_step1[1], parse_step1[2]);
}

std::ostream& operator<<(std::ostream& os, parsedModelTerm& p)
{
    os << p.shortModelTerm << ": funcName" << "[" << p.funcName << "] ";
    os << "variableString" << "[" << p.variableString << "] ";
    os << "variance" << "[" << p.allOptions["V"].valstring << "] ";
    os << "varianceStruct" << "[" << p.varianceStruct << "] ";
    os << "prior" << "[" << p.allOptions["prior"].valstring << "] ";
    os << "\n";
    return os;
}
