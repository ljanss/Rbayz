//
//  BayzR -- parVector.cpp
//

#include "parVector.h"
using Rsize_t = long int;

// common things for all contructors, this one is called at the end of every
// constructor because nelem must be set.
void parVector::common_constructor_items(parsedModelTerm & modeldescr, std::string namePrefix) {
   variables=modeldescr.variableString;                // as original, e.g A|B:C
   std::string tempname = modeldescr.variableString;   // and make parameter name from the variableString...
   size_t pos=0;                                       // that removes link-ID and changes :| to dots
   if( (pos=tempname.find('/')) != std::string::npos )
      tempname.erase(0, pos+1);
   if(namePrefix=="")
      Name = tempname;
   else
      Name=namePrefix + "." + tempname;
   while( (pos=Name.find_first_of(":|/",pos)) != std::string::npos) {
      Name[pos]='.';
      pos++;  // start re-search after currently replaced character
   }
   // [ToDo] Move disambiguation of parameter Name here??
   modelFunction=modeldescr.funcName;
   varianceStruct="-";
   val=Values.data;
   postMean.initWith(nelem,0.0l);
   postVar.initWith(nelem,0.0l);
   sumSqDiff.initWith(nelem, 0.0l);
   // [ToDo] parVector could also check the modeldescr for user-set 'traced' option
}

// contructor for par-vector with single element where the "variableString" is also used for the label
parVector::parVector(parsedModelTerm & modeldescr, double initval)
      : Values(), postMean(), postVar(), sumSqDiff() {
   nelem=1;
   Values.initWith(1, initval);
   Labels.push_back(modeldescr.variableString);
   common_constructor_items(modeldescr, "");
}

// constructor for single parameter value with prefix, used a.o. to make "var."
parVector::parVector(parsedModelTerm & modeldescr, double initval, std::string namePrefix)
      : Values(), postMean(), postVar(), sumSqDiff() {
   nelem=1;
   Values.initWith(1, initval);
   Labels.push_back(namePrefix + "." + modeldescr.variableString);
   size_t pos;  // still not good, now prefix falls out again ...
   if( (pos=Labels[0].find('/')) != std::string::npos )  // also done in param Name, so a bit
      Labels[0].erase(0, pos+1);                         // of same code in two places ...
   common_constructor_items(modeldescr, namePrefix);
}

// response model needs a constructor with a vector of values and vector of labels, and also uses namePrefix
parVector::parVector(parsedModelTerm & modeldescr, double initval, Rcpp::CharacterVector& inplabels,
            std::string namePrefix) : Values(), postMean(), postVar(), sumSqDiff() {
   nelem = inplabels.size();
   Values.initWith(nelem, initval);
   Labels.resize(inplabels.size());
   for(Rsize_t i=0; i<inplabels.size(); i++)
      Labels[i]=inplabels[i];
   common_constructor_items(modeldescr, namePrefix);
}

// many other model objects can initialize from a single scalar value and labels, the size
// needed is taken from labels size.
parVector::parVector(parsedModelTerm & modeldescr, double initval, Rcpp::CharacterVector& inplabels)
          : Values(), postMean(), postVar(), sumSqDiff() {
   nelem = inplabels.size();
   Values.initWith(nelem, initval);
   Labels.resize(inplabels.size());
   for(Rsize_t i=0; i<inplabels.size(); i++)
      Labels[i]=inplabels[i];
   common_constructor_items(modeldescr,"");
}

// nearly the same but labels is a vector<string>
parVector::parVector(parsedModelTerm & modeldescr, double initval, std::vector<std::string>& inplabels)
  : Values(), postMean(), postVar(), sumSqDiff() {
   nelem = inplabels.size();
   Values.initWith(nelem, initval);
   Labels.resize(inplabels.size());
   for(size_t i=0; i<inplabels.size(); i++)
      Labels[i]=inplabels[i];
   common_constructor_items(modeldescr,"");
}

// Update cumulative means and variances
void parVector::collectStats() {
   double olddev, newdev;
   this->count_collect_stats++;
   double n = double(count_collect_stats);
   if (count_collect_stats==1) {                  // at first sample collection store mean
      for(size_t i=0; i<nelem; i++) {
         this->postMean.data[i] = this->Values[i];
      }
   }
   else {                                     // can update mean, sumSqDiff and compute var
      for(size_t i=0; i<nelem; i++) {
         olddev = this->Values[i] - this->postMean.data[i]; // deviation with old mean
         this->postMean.data[i] += olddev/n;
         newdev = this->Values[i] - this->postMean.data[i]; // deviation with updated mean
         this->sumSqDiff.data[i] += olddev*newdev;
         this->postVar.data[i] = this->sumSqDiff.data[i]/(n-1.0l);
      }
   }
}

// Function to write (part of) parVector (for debugging purposes) - Rcout accepts this fine.
// operator<< overload must be defined as a non-member ...
std::ostream& operator<<(std::ostream& os, const parVector& p)
{
    os << p.Name << "[" << p.nelem << "] ";
    size_t loopsize = (p.nelem<5)? p.nelem : 5;
    for(size_t i=0; i<loopsize; i++) os << p.val[i] << " ";
    if(p.nelem >5 ) os << "...";
    os << "\n";
    return os;
}
