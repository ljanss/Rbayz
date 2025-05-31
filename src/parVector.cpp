//
//  BayzR -- parVector.cpp
//

#include "Rbayz.h"
#include "rbayzExceptions.h"
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
   for(pos=0; pos < tempname.size(); pos++) {
      if(tempname[pos]==':' || tempname[pos]=='|' || tempname[pos]=='/') tempname[pos]='.';
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
   // check trace option from the model-description
   // [ToDo]? This does not yet allow to switch off tracing where it is default on, to handle that,
   // need to check if the default is set before or after this parVector constructor ... switching it
   // off here where when the model constructor switches it back on as default does not work ...
   optionSpec trace_opt = modeldescr.allOptions["trace"];
   if (trace_opt.isgiven && trace_opt.valbool==true) {
      traced = 1;
      if(nelem>100) Rbayz::Messages.push_back("WARNING using 'trace' on "+Name+" (size="+std::to_string(nelem)+
             ") may need large memory; you could use 'save' instead to store samples in a file");
   }
   // check save option and open samples file if requested
   optionSpec save_opt = modeldescr.allOptions["save"];
   if(save_opt.isgiven && save_opt.valbool==true) {
      if( (openSamplesFile()) > 0) {
         throw generalRbayzError("Unable to open file for writing samples for " + Name);
      }
      saveSamples=true;
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
   count_collect_stats++;
   double n = double(count_collect_stats);
   if (count_collect_stats==1) {                  // at first sample collection store mean
      for(size_t i=0; i<nelem; i++) {
         postMean.data[i] = val[i];
      }
   }
   else {                                     // can update mean, sumSqDiff and compute var
      for(size_t i=0; i<nelem; i++) {
         olddev = val[i] - postMean.data[i]; // deviation with old mean
         postMean.data[i] += olddev/n;
         newdev = val[i] - postMean.data[i]; // deviation with updated mean
         sumSqDiff.data[i] += olddev*newdev;
         postVar.data[i] = sumSqDiff.data[i]/(n-1.0l);
      }
   }
}

int parVector::openSamplesFile() {
   std::string filename = "samples." + Name + ".txt";
   samplesFile = fopen(filename.c_str(),"w"); 
   return (samplesFile==0) ? 1 : 0; 
}

void parVector::writeSamples(int cycle) {
   // There are cases where samplesFile is not opened, it happens when saving with
   // other than the standard 'save' option (e.g. alpha's from the rn() model with kernel).
   if(saveSamples && samplesFile==0) {
      if ( (openSamplesFile()) > 0) {
         throw generalRbayzError("Unable to open file for writing samples for " + Name);
      }
   }
   if(saveSamples) {
      fprintf(samplesFile,"%d",cycle);
      for(size_t k=0; k<nelem; k++)
         fprintf(samplesFile," %g",Values[k]);
      fprintf(samplesFile,"\n");
   }
}

parVector::~parVector() {
   if(samplesFile != 0) fclose(samplesFile);
}

// Function to write name, size and first elements of a parVector (for debugging purposes),
// can be used with Rcpp::Rcout << `parVector`;
// ... operator<< overload must be defined as a non-member ...
std::ostream& operator<<(std::ostream& os, const parVector& p)
{
    os << p.Name << "[" << p.nelem << "] ";
    size_t loopsize = (p.nelem<5)? p.nelem : 5;
    for(size_t i=0; i<loopsize; i++) os << p.val[i] << " ";
    if(p.nelem >5 ) os << "...";
    return os;
}
