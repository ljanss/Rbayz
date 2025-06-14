#include <vector>
#include <string>
#include <Rcpp.h>
#include <cmath>
#include "Rbayz.h"
#include "parseFunctions.h"
#include "parsedModelTerm.h"
#include "modelBase.h"
#include "modelResp.h"
#include "modelMean.h"
#include "modelFixf.h"
#include "model_rn_ind.h"
#include "modelFreg.h"
#include "modelRreg.h"
#include "model_rn_cor.h"
#include "modelMixt.h"
#include "rbayzExceptions.h"
#include "simpleMatrix.h"
#include "simpleVector.h"
#include "modelVar.h"
#include "indepVarStr.h"
#include <unistd.h>

// [[Rcpp::plugins("cpp11")]]

// Definition of variable in Rbayz namespace; note these variables appear to persist between R calls,
// they are cleared / set again at the start of main().
std::vector<parVector**> Rbayz::parList;
std::vector<std::string> Rbayz::Messages;
bool Rbayz::needStop=false;
Rcpp::DataFrame Rbayz::mainData;

// [[Rcpp::export]]
Rcpp::List rbayz_cpp(Rcpp::Formula modelFormula, SEXP VE, Rcpp::DataFrame inputData,
                     Rcpp::IntegerVector chain, SEXP methodArg, int verbose,
                     Rcpp::Nullable<Rcpp::List> initVals_ = R_NilValue
                     )
//                   note VE and method are strings, it will be converted below
//                   note2: compiler also wants default values for any arguments placed after the
//                          Nullable initVals list, therefore the Nullable<list> comes last ...
{

   // clearing and setting global variables
   Rbayz::parList.clear();
   Rbayz::Messages.clear();
   Rbayz::needStop=false;
   Rbayz::mainData=inputData;

   // rbayz retains a small string describing last executed code that is sometimes added in errors
   std::string lastDone;

   // some variable are outside try{} because associated memory alloc needs to be cleaned up in catch{}
   modelResp* modelR = 0;
   std::vector<modelBase *> model;

   if (verbose > 0) Rcpp::Rcout << "R/bayz 0.10(.07)\n";

   try {     // normal execution builds a return list at the end of try{}; in case of
             // errors catch() builds a return list with the messages vector defined above.

      // split the modelFormula in a list of response (LHS) and explanatory (RHS) terms.
      std::string formulaAsCppstring = convertFormula(modelFormula);
      removeSpaces(formulaAsCppstring);
      std::vector<std::string> modelTerms         // response will be modelTerms[0] and forces in a
         = splitModelTerms(formulaAsCppstring);   // modelTerms[1] mn(0) or mn(1) for the intercept
      lastDone="Parsing model";
      if (verbose > 1) Rcpp::Rcout << "Parsing model done\n";

      // build response object - including response variance structure
      std::string VEstr =  Rcpp::as<std::string>(VE);
      parsedModelTerm parsedResponseVariable(modelTerms[0], VEstr);
      // here still need to add selecting different response objects based on variance structure
      modelR = new modelResp(parsedResponseVariable);   
      if (verbose > 1) Rcpp::Rcout << "Response model-object done\n";

      // Build vector of modelling objects from RHS terms (loop from term=1)
      if(verbose>2) Rcpp::Rcout << "Starting on building model objects ...\n";
      for(size_t term=1; term<modelTerms.size(); term++) {
         parsedModelTerm pmt(modelTerms[term]);
         if(verbose>2) Rcpp::Rcout << " ... building term " << term << " " << pmt.funcName << "()\n";
         if(pmt.funcName=="mn") {
            if(pmt.variableString=="1") model.push_back(new modelMean(pmt, modelR));
         }
         else if(pmt.funcName=="fx") {
            model.push_back(new modelFixf(pmt, modelR));
         }
         else if(pmt.funcName=="rn") {
            if(pmt.varianceStruct=="IDEN" || pmt.varianceStruct=="notgiven") {
               model.push_back(new model_rn_ind_iden(pmt, modelR));
            }
            else if (pmt.varianceStruct=="1kernel" || pmt.varianceStruct=="kernels") {
               model.push_back(new model_rn_cor_k0(pmt, modelR));
            }
            else {
               throw generalRbayzError("There is no class to model rn(...) with Variance structure " + pmt.allOptions["V"].valstring);
            }
         }
         else if (pmt.funcName=="rr") {
            // Note: varianceStruct DIAG, LASS, MIXT are always single (for now) - see parsedModelTerm.
            //       If there would be multiple with a DIAG it would be annotated as "mixed" varianceStruct.
            if(pmt.varianceStruct=="IDEN" || pmt.varianceStruct=="notgiven")
               model.push_back(new modelRregIden(pmt, modelR));
            else if (pmt.varianceStruct=="DIAG")
               model.push_back(new modelRregDiag(pmt, modelR));
            else if (pmt.varianceStruct=="LASS")
               model.push_back(new modelRregGRL(pmt, modelR));
            else if (pmt.varianceStruct=="MIXT") {
               modelRreg* rrmodel = new modelRregMixt(pmt, modelR);
               model.push_back(rrmodel);
               model.push_back(new modelMixt(pmt, rrmodel));
            }
            else
               throw generalRbayzError("There is no class to model rr(...) with Variance structure " + pmt.allOptions["V"].valstring);
         }
         else if (pmt.funcName=="rg") {    // [ToDo] work on adding rg() versions
            if(pmt.variablePattern=="onevar")
               model.push_back(new modelFreg(pmt, modelR));
//            else if (pmt.variablePattern="nestedreg")
//             need a new model object for the nested regression
            else
               throw generalRbayzError("Regression with the variable syntax " + pmt.variableString + " not yet supported\n");
         }
         else {
           throw generalRbayzError("Unknown model-function \'" + pmt.funcName + "\' at "+pmt.shortModelTerm);
         }
      } // end for(term ...) to build model
      lastDone="Model building";
      if (verbose>1) Rcpp::Rcout << "Model building done\n";
      if(Rbayz::needStop)
         throw(generalRbayzError("Quitting after model building because of errors"));

      // compute number of residuals (data points) and number of parameters in the model.
      // Response object is built first and parList[0] has residuals/fitted values.
      size_t nResiduals = (*(Rbayz::parList[0]))->nelem;
      size_t nParameters = 0;
      size_t nNAs = sum(modelR->missing);
      for(size_t i=1; i<Rbayz::parList.size(); i++) nParameters += (*(Rbayz::parList[i]))->nelem;
      {
         std::string s1="Note: data included total="+std::to_string(nResiduals)+" observed="+
                       std::to_string(nResiduals-nNAs)+" missing="+std::to_string(nNAs);
         std::string s2="Note: model build with "+std::to_string(nParameters)+" parameters";
         Rbayz::Messages.push_back(s1);
         Rbayz::Messages.push_back(s2);
      }
      if(verbose > 2) {
         Rcpp::Rcout << "Model-object overview (#, Name, Size, Traced, first Labels) after model building:\n";
         for(size_t i=0; i<Rbayz::parList.size(); i++) {
            Rcpp::Rcout << i << " " << (*(Rbayz::parList[i]))->Name << " " << (*(Rbayz::parList[i]))->nelem << 
                    " " << (*(Rbayz::parList[i]))->traced;
            Rcpp::Rcout << " " << (*(Rbayz::parList[i]))->Labels[0];
            if((*(Rbayz::parList[i]))->nelem > 1) Rcpp::Rcout << " " << (*(Rbayz::parList[i]))->Labels[1];
            if((*(Rbayz::parList[i]))->nelem > 2) Rcpp::Rcout << " ...";
            Rcpp::Rcout << "\n";
         }
      }

      // Make parameter name disambiguation - this skips residuals (parList[0])
      // I considered this can move to parVector where names are set.
/*
      {
         Rcpp::CharacterVector parNames;
         for(size_t i=1; i<Rbayz::parList.size(); i++)
            parNames.push_back((*Rbayz::parList[i])->Name);
         Rcpp::IntegerVector name_matches = Rcpp::match(parNames, parNames);
         std::vector<size_t> numb_matches(name_matches.size(),0);
         bool found_duplicates = FALSE;
         for(int i=0; i < name_matches.size(); i++) {
            if(name_matches[i] != (i+1))  {  // the i'th name matches (name_matches[i]-1)'th name in the list
               if(numb_matches[name_matches[i]-1] == 0) {
                  parNames[name_matches[i]-1] += ".1";
                  numb_matches[name_matches[i]-1]++; 
               }
               numb_matches[name_matches[i]-1]++;
               parNames[i] += "." + std::to_string(numb_matches[name_matches[i]-1]);
            }
            found_duplicates=TRUE;
         }
         if(found_duplicates) {  // copy new names back in parList->Names
            for(int i=0; i<parNames.size(); i++)
               (*(Rbayz::parList[i+1]))->Name = parNames[i];
         }
      }
      if (verbose>4) Rcpp::Rcout << "Parameter names disambiguation checked\n";
*/

      // Load initial values if given. An easy start is to only allow init-values from a run
      // with the same model - so that parameter-names and sizes all align.
      if(initVals_.isNotNull()) {
         Rcpp::List initVals(initVals_);
         Rcpp::DataFrame old_parameters = Rcpp::as<Rcpp::DataFrame>(initVals["Parameters"]);
         Rcpp::CharacterVector old_par_names = old_parameters["Param"];
         Rcpp::IntegerVector old_par_sizes = old_parameters["Size"];
         // check alignment of names and sizes between old and current build model
         bool match=true;
         if((size_t) old_par_sizes.size() == Rbayz::parList.size()){
            for(size_t i=0; i<Rbayz::parList.size(); i++) {
               if( ! (old_par_names[i] == (*(Rbayz::parList[i]))->Name && 
                        (size_t) old_par_sizes[i] == (*(Rbayz::parList[i]))->nelem) ) match=false;
            }
         }
         else
            match=false;
         if(match) {
            Rcpp::List old_estimates = initVals["Estimates"];
            for(size_t par=0; par<Rbayz::parList.size(); par++) {                     // for now loading fitted values from Estimates list,
               Rcpp::DataFrame par_data = Rcpp::as<Rcpp::DataFrame>(old_estimates[par]);    // but they are also stored in Residuals. 
               Rcpp::NumericVector par_pm = par_data["PostMean"];
               for(size_t row=0; row < (*(Rbayz::parList[par]))->nelem; row++)
                  (*(Rbayz::parList[par]))->val[row] = par_pm[row];
            }
            modelR->readjResid();  // residuals need to be reset to match loaded fitted values.
            modelR->restart();
            for(size_t mod=0; mod<model.size(); mod++) model[mod]->restart();
         }
         else {  // no match
            throw (generalRbayzError("Initialisation values cannot be used because names or sizes don't match"));
         }      // this could also be a warning, but there is no nice way to count and handle warnings
         Rcpp::Rcout << "Chain has been initialized with previous estimates\n";
      }
      if (verbose>4) Rcpp::Rcout << "Passed check init values\n";

      // Check the chain settings and find number of output samples by making list of output cycle-numbers.
      if (chain[0]==0 && chain[1]==0 && chain[2]==0) {  // chain was not set
         chain[0]=1100; chain[1]=100; chain[2]=10;
         Rcpp::Rcout << "Warning: chain was not set, running 1100 cycles but it may be too short for many analyses\n";
      }
      if (chain.size() != 3) throw (generalRbayzError("The chain settings do not have 3 elements"));
      if (chain[0] <= 0) throw (generalRbayzError("The chain length is zero or negative"));
      if (chain[1] < 0 || chain[2]<0 ) throw (generalRbayzError("The chain burnin or skip is negative"));
      if (chain[2] == 0) chain[2]=1;  // interpret skip zero as outputting all cycles
      Rcpp::IntegerVector outputCycleNumbers;
      for (int cycle=1; cycle <= chain[0]; cycle++) {         // exactly the same loop and condition to make
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0))  // output samples below
            outputCycleNumbers.push_back(cycle);
      }
      size_t nSamples = outputCycleNumbers.size();
      if (nSamples==0) throw (generalRbayzError("The chain settings do not make any output"));
      if (verbose>4) Rcpp::Rcout << "Chain checks done\n";

      // Find the number of traced parameters and set-up matrix to store samples of traced parameters
      size_t nTracedParam=0;
      if(verbose>0) Rcpp::Rcout << "Saving full traces for:";
      for(size_t i=0; i<Rbayz::parList.size(); i++) {
         if( (*(Rbayz::parList[i]))->traced ) {
            nTracedParam += (*(Rbayz::parList[i]))->nelem;
            if(verbose>0) {
               if((*(Rbayz::parList[i]))->nelem == 1) Rcpp::Rcout << " " << (*(Rbayz::parList[i]))->Labels[0];
               else {
                  for(size_t j=0; j< (*(Rbayz::parList[i]))->nelem; j++)
                     Rcpp::Rcout << " " << (*(Rbayz::parList[i]))->Name << "." << (*(Rbayz::parList[i]))->Labels[j];
               }
            }
         }
      }
      if(verbose>0) Rcpp::Rcout << "\n";
      Rcpp::NumericMatrix tracedSamples(nSamples,nTracedParam);

      // to show convergence set "nShow" interval and make vector to hold previously shown solutions
      int nShow = chain[0]/10;
      if( nShow < 1) nShow=1;
      int collect_first_conv = chain[0]/20;
      if (collect_first_conv < 1) collect_first_conv=1;
      Rcpp::NumericVector prevShowConv(int(nTracedParam), 0.5l);
      lastDone="Preparing to run MCMC";
      if (verbose>1) Rcpp::Rcout << "Preparing to run MCMC done\n";

      // Run the MCMC chain for method "Bayes" and "BLUPMC"
      // ------------------

      std::string method = Rcpp::as<std::string>(methodArg);
      if (method=="Bayes" || method=="BLUPMC") {
         if(verbose>0) {
            Rcpp::Rcout << "Cycle avgChange\n";
         }
         for (int cycle=1, save=0; cycle <= chain[0]; cycle++) {
            modelR->sample();
            for(size_t mt=0; mt<model.size(); mt++) model[mt]->sample();
            if(method=="Bayes") {
               modelR->sampleHpars();
               for(size_t mt=0; mt<model.size(); mt++) model[mt]->sampleHpars();
            }
            // At the 'skip' intervals and after burn-in:
            // 1) update posterior statistics using collectStats();
            // 2) save MCMC samples in memory for the 'traced' parameters;
            // 3) save MCMC samples on disk for parameters with 'saveSamples' option
            if ( (cycle > chain[1]) && (cycle % chain[2] == 0) ) {
               for(size_t mt=0; mt<model.size(); mt++) model[mt]->prepForOutput();
               for(size_t i=0; i<Rbayz::parList.size(); i++) (*(Rbayz::parList[i]))->collectStats();
               for(size_t i=0, col=0; i<Rbayz::parList.size(); i++) {
                  if( (*(Rbayz::parList[i]))->traced ) {
                     for(size_t j=0; j< (*(Rbayz::parList[i]))->nelem; j++) {
                        tracedSamples(save,col) = (*(Rbayz::parList[i]))->val[j]; col++;
                     }
                  }
                  if( (*(Rbayz::parList[i]))->saveSamples ) {
                     (*(Rbayz::parList[i]))->writeSamples(cycle);
                  }
               }
               save++;  // save is counter for output (saved) cycles
            }
            // at 'collect_first_conv' cycle store parameter values in prevShowConv to allow computing
            // first convergence; then at 'nShow' intervals show convergence on screen (when verbose > 0).
            if (cycle == collect_first_conv) {
               for(size_t i=0, col=0; i<Rbayz::parList.size(); i++) {
                  if( (*(Rbayz::parList[i]))->traced ) {
                     for(size_t j=0; j< (*(Rbayz::parList[i]))->nelem; j++) {                        
                        prevShowConv[col] = (*(Rbayz::parList[i]))->val[j];
                        col++;
                     }
                  }
               }
            }
            if (cycle % nShow == 0 && verbose>0) {
               Rcpp::Rcout << cycle;
               double conv_change=0.0l, conv_denom=0.0l, postmean;
               for(size_t i=0, col=0; i<Rbayz::parList.size(); i++) {
                  if( (*(Rbayz::parList[i]))->traced ) {
                     for(size_t j=0; j< (*(Rbayz::parList[i]))->nelem; j++) {
                        if (save==0) postmean = (*(Rbayz::parList[i]))->val[j];  // if nothing saved yet, using sampled
                        else postmean = (*(Rbayz::parList[i]))->postMean[j];     // value instead of postmean.
                        conv_change += abs(prevShowConv[col] - postmean);
                        conv_denom += abs(prevShowConv[col]);
                        prevShowConv[col] = postmean;
                        col++;
                     }
                  }
               }
               conv_change /= conv_denom;
               Rcpp::Rcout << " " << conv_change << "\n";
            }
         } // end for(cycle ...)
      }

/*    else if (method=="BLUP") {     // insert here BLUP version
      }
*/

      lastDone="Finished running MCMC";
      if (verbose>1) Rcpp::Rcout << "Finished running MCMC\n";

      // Build tables to go in the output:
      // ---------------------------------

      // 1. "Parameter" information table
      Rcpp::CharacterVector parNames, parModelFunc, parVariables, parVarStruct;
      Rcpp::IntegerVector parSizes, parTraced;
      for(size_t i=0; i<Rbayz::parList.size(); i++) {
         parNames.push_back((*(Rbayz::parList[i]))->Name);
         parModelFunc.push_back((*(Rbayz::parList[i]))->modelFunction);
         parVariables.push_back((*(Rbayz::parList[i]))->variables);
         parVarStruct.push_back((*(Rbayz::parList[i]))->varianceStruct);
         parSizes.push_back((*(Rbayz::parList[i]))->nelem);
         parTraced.push_back((*(Rbayz::parList[i]))->traced);
      }
      Rcpp::DataFrame parInfo = Rcpp::DataFrame::create
               (Rcpp::Named("ModelTerm")=parModelFunc, Rcpp::Named("Variables")=parVariables, 
                Rcpp::Named("Param")=parNames,Rcpp::Named("Variance")=parVarStruct,
                Rcpp::Named("Size")=parSizes, Rcpp::Named("Traced")=parTraced);
//      parInfo.attr("row.names") = parNames;

      // 2. "Estimates": now a list with a data frame for each parameter-vector
      Rcpp::List estimates = Rcpp::List::create();
      for(size_t i=0; i < Rbayz::parList.size(); i++) {    // For the moment including fitval from parList[0], because init
         size_t nr = (*(Rbayz::parList[i]))->nelem;        // values reads it from there, but they are also stored in "Residuals" ...
         Rcpp::CharacterVector rcpp_labels(nr);
         Rcpp::NumericVector rcpp_postmeans(nr);
         Rcpp::NumericVector rcpp_postSDs(nr);
         for(size_t row=0; row<nr; row++) {
            rcpp_labels[row] = (*(Rbayz::parList[i]))->Labels[row];
            rcpp_postmeans[row] = (*(Rbayz::parList[i]))->postMean[row];
            rcpp_postSDs[row] = sqrt( (*(Rbayz::parList[i]))->postVar[row] );
         }
         Rcpp::DataFrame thispar_estimates = Rcpp::DataFrame::create
              (Rcpp::Named("Level")=rcpp_labels, Rcpp::Named("PostMean")=rcpp_postmeans,
              Rcpp::Named("PostSD")=rcpp_postSDs);
         estimates.push_back(thispar_estimates,(*(Rbayz::parList[i]))->Name);
      }
      lastDone="Computing postMeans and PostSDs";

      // 3. "Samples" table (this is the matrix tracedSamples with row and col-names added)
      Rcpp::CharacterVector sampleRowNames = Rcpp::wrap(modelR->par->Labels);
      Rcpp::CharacterVector sampleColNames;
      for(size_t i=0; i<Rbayz::parList.size(); i++) {
         if( (*(Rbayz::parList[i]))->traced ) {
            if( (*(Rbayz::parList[i]))->nelem==1) sampleColNames.push_back( (*(Rbayz::parList[i]))->Name);
            else {
               std::string s = (*(Rbayz::parList[i]))->Name;
               for(size_t j=0; j< (*(Rbayz::parList[i]))->nelem; j++)
                  sampleColNames.push_back(s + (*(Rbayz::parList[i]))->Labels[j]);
            }
         }
      }
      Rcpp::colnames(tracedSamples) = sampleColNames;
      Rcpp::rownames(tracedSamples) = Rcpp::as<Rcpp::CharacterVector>(outputCycleNumbers); 
      /* I couldn't get this colnames() and rownames() working, it gives a compiler error that the Rcpp::NumericMatrix
         can't be conveted to SEXP object - but online examples show this should work ...
         */
      // work-around to set row and colnames of a matrix: add "dimnames" attribute as a list of 2
//      Rcpp::List tracedSamplesNames = Rcpp::List::create(sampleRowNames,sampleColNames);
//      tracedSamples.attr("dimnames") = tracedSamplesNames;

      // 4. "Residuals" table: compute residuals from Y.data and fitted value (par-vector in modelR).
      // The 'resid' in modelR cannot be used because that one is a sampled state, not a posterior mean.
      // Need to think about modifications for non-linear models, then residual may also need to be stored
      // and averaged because it is more difficult to computer from the Y and fitted value?
      Rcpp::NumericVector fitval(nResiduals);
      Rcpp::NumericVector resid(nResiduals);
      for(size_t i=0; i<nResiduals; i++) {
         fitval[i] = modelR->par->postMean[i];
         if(modelR->missing[i]) resid[i] = NA_REAL;  // residual NA for missing data, but fitval exists!
         else resid[i] = modelR->Y.data[i] - modelR->par->postMean[i];
      }
      Rcpp::NumericVector residuals(resid);
      Rcpp::CharacterVector residRowNames = Rcpp::wrap(modelR->par->Labels);
      residuals.names() = residRowNames;

      // Build the final return list
      Rcpp::List result = Rcpp::List::create();
      result.push_back(0,"nError");
      if(Rbayz::Messages.size()>0) 
         result.push_back(Rbayz::Messages,"Messages");
      result.push_back(parInfo,"Parameters");
      result.push_back(tracedSamples,"Samples");
      result.push_back(estimates,"Estimates");
      result.push_back(residuals,"Residuals");
      result.push_back(chain,"Chain");
      lastDone="Filling return list";
      if (verbose>1) Rcpp::Rcout << "Ready filling return list\n";

      // clean-up and normal termination
      // ------------------
      if(modelR != 0) delete modelR;
      for(size_t i=0; i<model.size(); i++) delete model[i];
      return(result);

   } // end try{}

   // Catching (all kinds of) errors (and build return list with the error message)
   // ------------------------------

   catch (generalRbayzError &err) {
      Rbayz::Messages.push_back(err.what());
   }
   catch (std::exception &err) {
      std::string s = std::string(err.what()) + " after " + lastDone;
      Rbayz::Messages.push_back(s);
   }
   catch (...) {
      Rbayz::Messages.push_back("An unknown error occured in bayz after: "+lastDone);
   }
   // Build a return list that only has the error messages list.
   Rcpp::List result = Rcpp::List::create();
   result.push_back(Rbayz::Messages.size(),"nError");
   result.push_back(Rbayz::Messages,"Messages");
   Rcpp::Rcout << "Bayz finished with (last) error: " << Rbayz::Messages[Rbayz::Messages.size()-1] << std::endl;
   Rcpp::Rcout << "There may be more messages or errors - use summary() or check <output>$Errors to see all" << std::endl;
//   if(modelR != 0) delete modelR;
//   for(size_t i=0; i<model.size(); i++) delete model[i];
   return(result);

}  // end rbayz_cpp main function

