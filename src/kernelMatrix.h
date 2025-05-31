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

   kernelMatrix(varianceSpec var_descr)
         : labeledMatrix(), weights() {
//    note: kernelMatrix starts with an empty labeledMatrix, parent constructors have not done anything,
//    accept for having the matrix and vectors for storing data and row and column labels.
//    The initWith() used at the end to copy eigenvec contents in the object is not simpleMatrix' initWith,
//    but labeledMatrix' version that also handles copying row and column labels.

      Rcpp::NumericMatrix kerneldata = Rcpp::as<Rcpp::NumericMatrix>(var_descr.kernObject);
   	Rcpp::Function eig("eigen");
	   Rcpp::List eigdecomp;
   	try {
	   	eigdecomp = eig(kerneldata);
   	}
   	catch(std::exception &err) {
         throw(generalRbayzError("An error occurred running eigen(): "+std::string(err.what())));
   	}

   	Rcpp::NumericVector eigvalues = eigdecomp["values"];
	   Rcpp::NumericMatrix eigvectors = eigdecomp["vectors"];

      // re-attach the dimnames again to the eigvectors matrix for correct further processing
      if (kerneldata.hasAttribute("dimnames")) {
         Rcpp::List dimnames = Rcpp::as<Rcpp::List>(kerneldata.attr("dimnames"));
         eigvectors.attr("dimnames") = dimnames;
      }

      // Get / check / set dim_size (dim) and/or dim_pct (dimp) options
      double dim_pct=0;
      int dim_size=0;
      optionSpec dim_opt = var_descr["dim"];
      optionSpec dimp_opt = var_descr["dimp"];
      if( dim_opt.isgiven ) {
         dim_size = (int) dim_opt.valnumb[0];
         if(dim_size <= 0 || dim_size > eigvalues.size()) {
            Rbayz::Messages.push_back("Warning: invalid dim setting <" + std::to_string(dim_size) + "> processing kernel " + var_descr.keyw + ", setting default dimp=90");
            dim_size=0;  // if dim not well set this does not trigger error,
            dim_pct=90;  // but goes back to cutting off on 90% of variance.
         }
      }
      else {  // without 'dim' option, check for 'dimp' (note: dim will be used when both are set!)
         if( dimp_opt.isgiven ) {
            dim_pct = dimp_opt.valnumb[0];
            if(dim_pct <= 0 || dim_pct > 100) {
               Rbayz::Messages.push_back("Warning: invalid dimp setting <" + std::to_string(dim_pct) + "> processing kernel " + var_descr.keyw + ", setting default dimp=90");
               dim_pct=90;  // also here not error, but fall back to default dimp=90
            }
         }
         else  // no options set: take default
            dim_pct=90;
      }
      double sumeval = 0.0l;             // the sum of all positive eigenvalues and their count
      size_t counted_positive_evals=0;
      for (size_t i = 0; i < unsigned(eigvalues.size()) && eigvalues[i] > 0; i++) {
         sumeval += eigvalues[i];
         counted_positive_evals++;
      }
      if(dim_size==0) {        // need to get a dim_size from dim_pct
         double eval_cutoff = dim_pct * sumeval / 100.0l;
         double sum_part = 0.0l;
         while (sum_part < eval_cutoff) sum_part += eigvalues[dim_size++];
      }
      else {                 // dim_size is set, check var explained and that it does not cover negative evals
         if (dim_size > counted_positive_evals) dim_size = counted_positive_evals;
         double sum_part = 0.0l;
         for (size_t i = 0; i < dim_size; i++) {
            sum_part += eigvalues[i];
         }
         dim_pct = 100.0 * sum_part / sumeval;
      }
      std::string s = "Note: in " + var_descr.optionText + " for kernel " + var_descr.keyw + " using dimp=" + std::to_string(dim_pct)
                     + " and dim=" + std::to_string(dim_size);
      Rbayz::Messages.push_back(s);
      this->initWith(eigvectors, var_descr.keyw, dim_size);
      weights.initWith(eigvalues, dim_size);
   }

   ~kernelMatrix() {
   }

   // Add a kernel (make the kronecker product) to the stored kernel in the object.
   void addKernel(kernelMatrix* K2) {
      // make list (and sort it) of all interaction evalues
 	   std::vector<double> evalint(this->ncol*K2->ncol, 0.0l);
      double sumeval = 0.0l;
	   for (size_t i = 0; i<this->ncol; i++) {
		  for (size_t j = 0; j<K2->ncol; j++) {
           evalint[i*K2->ncol + j] = this->weights[i] * K2->weights[j];
			  sumeval += evalint[i*K2->ncol + j];
		  }
      }
	   std::sort(evalint.begin(), evalint.end(), std::greater<double>());
      // Determine cut-off eval to reach rrankpct cumulative sum
      double rrankpct=90;          // rrankpct not coming correctly from the modelterm now
	   double eval_sum_cutoff = rrankpct * sumeval / 100.0l;  // this is cut-off on cumulative eval
	   sumeval = 0.0l;
	   unsigned long nEvalUsed=0;
	   while ((sumeval += evalint[nEvalUsed]) < eval_sum_cutoff) nEvalUsed++;
	   double eval_min_cutoff = evalint[nEvalUsed];  // this is cut-off evalue to get rrankpct cumulative
      nEvalUsed++;
      // Set-up a temporary  matrix and vector to make the evectors and evalues for the interaction kernel.
      // For the moment I keep setting up for all combinations. 
      size_t nLevel1 = this->nrow;
      size_t nLevel2 = K2->nrow;
      simpleMatrix tempEvecs(nLevel1*nLevel2, nEvalUsed);
      simpleDblVector tempEvals(nEvalUsed);
      size_t k = 0;
      for(size_t i=0; i<this->ncol; i++) {
         for(size_t j=0; j<K2->ncol; j++) {
            if ( this->weights[i]*K2->weights[j] >= eval_min_cutoff ) {
               tempEvals.data[k] = this->weights[i] * K2->weights[j];
               for(size_t rowi=0; rowi<nLevel1; rowi++) {
                  for(size_t rowj=0; rowj<nLevel2; rowj++) {
                     tempEvecs.data[k][i*nLevel2+j] = this->data[i][rowi] * K2->data[j][rowj];
                  }
               }
               k++;
            }
         }
      }
      // Make row-colnames for the new combined matrix
      std::vector<std::string> tempRownames;
      tempRownames.reserve(nLevel1*nLevel2);
      for(size_t rowi=0; rowi<nLevel1; rowi++) {
         for(size_t rowj=0; rowj<nLevel2; rowj++) {
            tempRownames.push_back(this->rownames[rowi]+"."+K2->rownames[rowj]);
         }
      }
      std::vector<std::string> tempColnames =  generateLabels("col",nEvalUsed);
      Rcpp::Rcout << "Interaction kernel retains " << nEvalUsed << " eigenvectors\n";
      // Swap the old data with the new data.
      // Note the old data is removed when tempEvecs, tempEvals and tempLabels here go out of scope.
      this->swap(&tempEvecs);
      weights.swap(&tempEvals);
      std::swap(rownames,tempRownames);
      std::swap(colnames,tempColnames);
   }

   // [ToDo] name weights is not so good, this is a class holding an eigendecomp, it is
   // quite ok to call it eigen values. Weights becomes confusing in other parts of the code.
   simpleDblVector weights;
   
};

#endif /* kernelMatrix_h */
