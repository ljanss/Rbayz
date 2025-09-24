//  kernelMatrix.h
//  Storage of kernel as eigen-decomposition.
//  Derives from simpleMatrix using the simpleMatrix() empty contructor, so that
//  no matrix data is stored yet from the parent constructor.
//
//  Created by Luc Janss on 06/05/2021.
//

#ifndef kernelMatrix_h
#define kernelMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "Rbayz.h"
#include "rbayzExceptions.h"
#include "parseFunctions.h"
#include <string>
#include <vector>
#include "labeledMatrix.h"
#include "simpleVector.h"
#include "nameTools.h"
#include "optionsInfo.h"

class kernelMatrix : public labeledMatrix {

public:

   kernelMatrix(varianceSpec var_descr); 
   kernelMatrix(varianceSpec var_descr, double dim_pct);
   ~kernelMatrix() {
   }

   // Add a kernel (make the kronecker product) to the stored kernel in the object.
   void addKernel(kernelMatrix* K2);

   // [ToDo] name weights is not so good, this is a class holding an eigendecomp, it is
   // quite ok to call it eigen values. Weights becomes confusing in other parts of the code.
   simpleDblVector weights;
   double sumEvalues;
   
};

#endif /* kernelMatrix_h */
