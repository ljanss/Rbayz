//
// BayzR -- dataCovar.h
// Storage class to hold a covariate.
// The constructors take an R-object, attempts to convert it to a numerical vector, and
// throw error when this fails. With an update, the original variable name needs to be passed
// too to allow printing a meaninful error message when conversion fails.
//
// Created by Luc Janss on 14/01/2020.
//

#ifndef dataCovar_h
#define dataCovar_h

#include <Rcpp.h>
#include "simpleVector.h"
#include "rbayzExceptions.h"

class dataCovar : public simpleDblVector {

public:

   // Default one centering and filling NA
   dataCovar(Rcpp::RObject x, std::string varname) : simpleDblVector() {
      try {
         Rcpp::NumericVector temp = Rcpp::as<Rcpp::NumericVector>(x);
         Rcpp::LogicalVector missing = Rcpp::is_na(temp);
         initWith(temp);
         size_t countnoNA=0;
         double sumnoNA=0.0l;
         for(size_t row=0; row < unsigned(temp.size()); row++) {
            if(!missing[row]) {
               sumnoNA += data[row];
               countnoNA++;
            }
         }
         offset = sumnoNA / double(countnoNA);
         for(size_t row=0; row < unsigned(temp.size()); row++) {
            if(missing[row])
               data[row] = 0.0;
            else
               data[row] -= offset;
         }
      }
   	catch(std::exception &err) {
         throw(generalRbayzError("Error in dataCovar transforming variable " + varname + 
                                 " to a numerical vector: "+std::string(err.what())));
   	}
   }

   // constructor with optional center and filling NA
   dataCovar(Rcpp::RObject x, bool center, bool fillNA, std::string varname) : simpleDblVector() {
      try {
         Rcpp::NumericVector temp = Rcpp::as<Rcpp::NumericVector>(x);
         Rcpp::LogicalVector missing = Rcpp::is_na(temp);
         initWith(temp);
         size_t countnoNA=0;
         double sumnoNA=0.0l;
         for(size_t row=0; row < unsigned(temp.size()); row++) {
            if(!missing[row]) {
               sumnoNA += data[row];
               countnoNA++;
            }
         }
         if(countnoNA < nelem && !fillNA)
            throw(generalRbayzError("Error in dataCovar: covariate " + varname + 
                                    " has missing values which is not allowed here"));
         double meancovar = sumnoNA / double(countnoNA);
         if (center) {
            offset = meancovar;   // store meancovar in 'offset' only when data is centered
            for(size_t row=0; row < unsigned(temp.size()); row++) {
               if(!missing[row]) data[row] -= meancovar;
            }
         }
         if (fillNA && center) {
            for(size_t row=0; row < unsigned(temp.size()); row++) {
               if(missing[row]) data[row] = 0;
            }
         }
         if (fillNA && !center) {
            for(size_t row=0; row < unsigned(temp.size()); row++) {
               if(missing[row]) data[row] = meancovar;
            }
         }
      }
   	catch(std::exception &err) {
         throw(generalRbayzError("Error in dataCovar transforming variable " + varname + 
                                 " to a numerical vector: "+std::string(err.what())));
   	}
   }

   ~dataCovar() {
   }

   // when fitting a regression on covar, there will be a part offset * regcoeff missing from
   // the model intercept because that part has been substracted in the centering of the covar.
   double offset=0.0l;
   
};

#endif /* dataCovar_h */
