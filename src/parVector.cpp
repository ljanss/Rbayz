//
//  BayzR -- parVector.cpp
//

#include "Rbayz.h"
#include "parVector.h"
using Rsize_t = long int;

// common things for all contructors, this one is called at the end of every
// constructor because nelem must be set.
void parVector::common_constructor_items(parsedModelTerm & modeldescr, std::string namePrefix) {
   variables=modeldescr.variableString;                // as original, e.g A|B:C
   std::string tempname = modeldescr.variableString;   // and make parameter name from the variableString...
   size_t pos=0;                                       // that removes link-ID and changes :| to dots
   if( (pos=tempname.find('/')) != std::string::npos ) tempname.erase(0, pos+1);
   if(namePrefix!="")
      tempname = namePrefix + "." + tempname;
   while( (pos=tempname.find_first_of(":|/",pos)) != std::string::npos) {
      tempname[pos]='.';
      pos++;  // start re-search after currently replaced character
   }
   Name=tempname;
   // Check for duplicate names in already stored parameter-list;
   // The last element in parList should be the one currently constructed, because
   // modelBase constructor runs first to put the par-pointer in parList, after which
   // the 'downstream' constructors run to allocate the parVector ... so run the loop
   // until < parList.size()-1. Also skip element 0 in parList (it is residuals).
   {
      bool duplicateFound;
      int duplicateCount=0;
      do {
         duplicateFound=false;
         for(size_t i=1; i < (Rbayz::parList.size()-1); i++) {
            if((*Rbayz::parList[i])->Name == tempname) {
               duplicateFound=true;
               break;
            }
         }
         if(duplicateFound) {
            duplicateCount++;
            tempname=Name+std::to_string(duplicateCount);
         }
      }
      while(duplicateFound);
   }
   Name=tempname;
   modelFunction=modeldescr.funcName;
   varianceStruct="-";
   val=Values.data;
   postMean.initWith(nelem,0.0l);
   postVar.initWith(nelem,0.0l);
   sumSqDiff.initWith(nelem, 0.0l);
   // get save and trace options from the model-description
   std::string traceopt = modeldescr.options["trace"];
   if(traceopt=="TRUE" || traceopt=="T" || traceopt=="1" || traceopt=="y") {
      traced = 1;
      if(nelem>100) Rbayz::Messages.push_back("WARNING using 'trace' on "+Name+" (size="+std::to_string(nelem)+
             ") may need large memory; you could use 'save' instead to store samples in a file");
   }
   else if (traceopt=="FALSE" || traceopt=="F" || traceopt=="0" || traceopt=="n") {
      traced = 0;
   }
   else {
      // error
   }
   std::string saveopt = modeldescr.options["save"];
   if(saveopt=="TRUE" || saveopt=="T" || saveopt=="1" || saveopt=="y") {
   }
   else if (saveopt=="FALSE" || saveopt=="F" || saveopt=="0" || saveopt=="n") {
   }
   else {
      // error
   }
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
   std::string templabel = modeldescr.variableString;
   size_t pos;
   if( (pos=templabel.find('/')) != std::string::npos ) templabel.erase(0, pos+1);   
   Labels.push_back(namePrefix + "." + templabel);
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

// parVector constructor that uses another 'related' parVector to copy size and labels, and name it as
// namePrefix + name of the relatedPar. This is used now in modelHelper class to add objects to hold
// additional parameter vectors.
parVector::parVector(parsedModelTerm & modeldescr, double initval, parVector & relatedPar, std::string namePrefix)
         : Values(), postMean(), postVar(), sumSqDiff() {
   nelem = relatedPar.nelem;
   Values.initWith(nelem, initval);
   Labels.resize(nelem);
   for(size_t i=0; i<nelem; i++)
      Labels[i]=relatedPar.Labels[i];
   common_constructor_items(modeldescr,namePrefix);
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
