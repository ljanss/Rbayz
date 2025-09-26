//  matrixClasses.cpp
//  Code for labeledMatrix and kernelMatrix classes

#include "labeledMatrix.h"
#include "rbayzExceptions.h"
#include "nameTools.h"

// ----------------- labeledMatrix class --------------------

// add/copy names from Rcpp matrix in the labeledMatrix object.
// Throws errror if rownames not available, auto-fills colnames if colnames not available
// This is now a member function of labeledMatrix, so that the object can call addRowNames
// on itself to get its row or colnames filled.
void labeledMatrix::addRowColNames(Rcpp::NumericMatrix M, std::string & name) {
   rownames = getMatrixNames(M, 1);
   if(rownames.size()==0) {  // rownames empty not allowed
      throw generalRbayzError("No rownames on matrix " + name + "\n");
   }
   colnames = getMatrixNames(M, 2);
   if (colnames.size()==0) { // colnames empty, fill auto colnames
      colnames = generateLabels("col",M.ncol());
   }
}

// Version restricting the column labels copied to 'useCol' columns, part of constructor
// with 'useCol' limit.
// Proper setting of useCol (>=1 and <= M.ncol) is not tested here, it is assumed that this
// labeling function is only used in labeledMatrix contructor, and then simpleMatrix() constructor
// will throw errors for improper setting of useCol.
void labeledMatrix::addRowColNames(Rcpp::NumericMatrix M, std::string & name, size_t useCol) {
   rownames = getMatrixNames(M, 1);
   if(rownames.size()==0) {  // rownames empty not allowed
      throw generalRbayzError("No rownames on matrix " + name + "\n");
   }
   std::vector<std::string> tempnames = getMatrixNames(M, 2);
   if(tempnames.size()==0) {   // no colnames supplied, auto-generate useCol colnames
      colnames = generateLabels("col",useCol);
   }
   else {                     // copy useCol columns from tempnames in colnames labels
      colnames.resize(useCol);
      for(size_t i=0; i<useCol; i++) colnames[i] = tempnames[i];
   }
}

labeledMatrix::labeledMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(col) {
   // Need to temporarily redo the conversion of the input Robject to
   // Rcpp::NumericMatrix to retrieve row and col names.
   Rcpp::NumericMatrix Rmatrix = Rcpp::as<Rcpp::NumericMatrix>(col);
   this->addRowColNames(Rmatrix, name);
}

void labeledMatrix::initWith(Rcpp::NumericMatrix & M, std::string & name, size_t useCol) {
   simpleMatrix::initWith(M, useCol);
   this->addRowColNames(M, name, useCol);
}

// ----------------- kernelMatrix class --------------------

// kernelMatrix constructor with no dim_pct setting, this calls the other constructor,
// setting the default dim_pct = 90
kernelMatrix::kernelMatrix(varianceSpec var_descr) : kernelMatrix(var_descr, 90.0) {
}

// kernelMatrix constructor. The dim_pct is used as a default; when there are options dim or dimp
// in the kernel specification, these are used instead of the default dim_pct.
kernelMatrix::kernelMatrix(varianceSpec var_descr, double dim_pct) : labeledMatrix(), weights() {

// note: kernelMatrix starts with an empty labeledMatrix, parent constructors have not done anything,
// accept for having the matrix and vectors for storing data and row and column labels.
// The initWith() used at the end to copy eigenvec contents in the object is not simpleMatrix' initWith,
// but labeledMatrix' version that also handles copying row and column labels.

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
   sumEvalues = 0.0l;             // the sum of all positive eigenvalues and their count
   size_t counted_positive_evals=0;
   for (size_t i = 0; i < unsigned(eigvalues.size()) && eigvalues[i] > 0; i++) {
      sumEvalues += eigvalues[i];
      counted_positive_evals++;
   }
   if(dim_size==0) {        // need to get a dim_size from dim_pct
      double eval_cutoff = dim_pct * sumEvalues / 100.0l;
      double sum_part = 0.0l;
      while (sum_part < eval_cutoff) sum_part += eigvalues[dim_size++];
   }
   else {                 // dim_size is set, check var explained and that it does not cover negative evals
      if (dim_size > counted_positive_evals) dim_size = counted_positive_evals;
      double sum_part = 0.0l;
      for (size_t i = 0; i < dim_size; i++) {
         sum_part += eigvalues[i];
      }
      dim_pct = 100.0 * sum_part / sumEvalues;
   }
   std::string s = "Note: in " + var_descr.optionText + " for kernel " + var_descr.keyw + " using dimp=" + std::to_string(dim_pct)
                  + " and dim=" + std::to_string(dim_size);
   Rbayz::Messages.push_back(s);
   this->initWith(eigvectors, var_descr.keyw, dim_size);
   weights.initWith(eigvalues, dim_size);
}

/*
'mergeKernels': make kronecker product of 'this' object kernel and a second one, and replace 'this'
with the merged one. This also updates the factor indexing and labels to match the interactions of the
two kernels.
I'd say the code should work to repeat this for multiple kernels, but I have not tested that yet.
An older version was making a selection of evecs to keep (see oldcode folder), but now it simply uses
all combinations! Therefore it should be used with caution, as the merged kernel can become very large;
the calling code (in model_rn_cor_k0) overviews how many kernels are merged and computes total size of
the final result, and checks if that fits in the maxmem allowed memory.
*/
void kernelMatrix::addKernel(kernelMatrix* K2) {
   size_t nLevel1 = this->nrow;
   size_t nLevel2 = K2->nrow;
   size_t nColumn1 = this->ncol;
   size_t nColumn2 = K2->ncol;
   simpleMatrix tempEvecs(nLevel1*nLevel2, nColumn1*nColumn2);
   simpleDblVector tempEvals(nColumn1*nColumn2);
   std::vector<std::string> tempColnames(nColumn1*nColumn2,"");
   size_t k = 0;
   for(size_t i=0; i<this->ncol; i++) {
      for(size_t j=0; j<K2->ncol; j++) {
         k = i*nColumn2 + j; // map column i,j to new column k
         tempEvals.data[k] = this->weights[i] * K2->weights[j];
         tempColnames[k]=this->colnames[i]+"."+K2->colnames[j];
         for(size_t rowi=0; rowi<nLevel1; rowi++) {
            for(size_t rowj=0; rowj<nLevel2; rowj++) {
               tempEvecs.data[k][i*nLevel2+j] = this->data[i][rowi] * K2->data[j][rowj];
            }
         }
      }
   }
   // Make rownames for the new combined matrix
   std::vector<std::string> tempRownames;
   tempRownames.reserve(nLevel1*nLevel2);
   for(size_t rowi=0; rowi<nLevel1; rowi++) {
      for(size_t rowj=0; rowj<nLevel2; rowj++) {
         tempRownames.push_back(this->rownames[rowi]+"."+K2->rownames[rowj]);
      }
   }
   // Swap the old 'this' data with the new data; the old data in 'this' kernel is then in the temp
   // objects and will be deleted when the temps go out of scope at the end of this function.
   // The data that was in K2 must be deleted by the calling code.
   this->swap(&tempEvecs);
   weights.swap(&tempEvals);
   std::swap(rownames,tempRownames);
   std::swap(colnames,tempColnames);
}

