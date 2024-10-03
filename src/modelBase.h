//  modelBase.h
//
//  Base class for model (computational) classes.
//  This defines common interface with par-vector and sample() vector, many
//  other details specific for response, explanatory variables or variances come
//  in derived classes.
//
//  Created by Luc Janss on 03/08/2018.

#ifndef modelBase_h
#define modelBase_h

#include <Rcpp.h>
#include <vector>
#include <stdio.h>
#include "parVector.h"
//#include <unistd.h>

extern std::vector<parVector**> parList;

class modelBase {
   
public:

   modelBase() {
      parList.push_back(&par);
   };

   virtual ~modelBase() {
   }
   
   virtual void sample() = 0; 

   // prepForOutput is for model classes that need to make a transform of
   // parameters for output, the base class defines an 'empty' version.
   virtual void prepForOutput() { };

   // saveSamples saves MCMC samples in a file, it will kick in when the constructor
   // has set save_samples=true. There is a definition in modelClasses.cpp file.
   virtual int openSamplesFile();
   virtual void saveSamples();
   bool saveSamples = false;
   FILE* samplesFile;

   parVector* par=0;

};

#endif /* modelBase_h */
