//  dataMatrix.h
//  Class for storage of matrix of covariates (e.g. input in rr() model).
//  Derives from labeledMatrix so it has rownames (and optionally colnames).
//  Added features compared to labeledMatrix: makes column centering, allows NA
//  in input data, which are replaced with column means; adds a member variable
//  sumColVar that stores the sum of column variances.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataMatrix_h
#define dataMatrix_h

#include <math.h>  // needed?
#include "labeledMatrix.h"

class dataMatrix : public labeledMatrix {

public:
   dataMatrix(Rcpp::RObject col, std::string & name) : labeledMatrix(col, name) {
      // column-center the matrix data and fill missings with column mean.
      // (and after centering the column means are zero for all columns).
      double * datacol;
      size_t i,j, nobs;
      double sum;
      sumColVar = 0.0l;
      for(i=0; i<ncol; i++) {
         datacol = data[i];
         sum=0.0l; nobs=0;
         double colVar = 0.0l;
         for(j=0; j<nrow; j++) {
            if (! std::isnan(datacol[j]) ) { // cannot catch it as NA_REAL after is has been copied
               sum+=datacol[j];         // in the simpleMatrix structure ... maybe can be improved?
               nobs++;
            }
         }
         sum /= double(nobs);
         for(j=0; j<nrow; j++) {
            if ( std::isnan(datacol[j]))
               datacol[j] = 0.0l;
            else 
               datacol[j] -= sum;
            colVar += datacol[j]*datacol[j];
         }
         sumColVar += colVar/double(nrow);
      }
   }

   ~dataMatrix() {
   }

   double sumColVar;

};

#endif /* dataMatrix_h */
