//
//  BayzR --- modelRreg.hpp
//
//  Computational class to model random regressions on a matrix of covariates.
//  This uses methods from modelMatrix, only some parameter names need to be set.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRreg_h
#define modelRreg_h

#include <Rcpp.h>
#include "modelMatrix.h"
#include "indepVarStr.h"
#include "dataMatrix.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelMatrix(modeldescr, rmod)   {
      if(checkOptions(modeldescr.options, "V prior pvals")>0) {
         throw(generalRbayzError("ERROR: unrecognized option(s) in "+modeldescr.shortModelTerm));
      }
   }

   ~modelRreg() {
      delete varmodel;
   }
   
   // sample() now updated to handle zero variance (inf weight) from the variance model:
   // - decorrect() runs to correct for old estimate, but skips if old one was already zero
   // - for inf weight, regcoeff is set to zero, else it will go through the sampling procedure
   // - correct() runs in case regcoeff != zero, but skips if it is zero
   void sample() {
      double inf = std::numeric_limits<double>::infinity();
      double lhs, rhs;
      for(size_t k=0; k < M->ncol; k++) {
         resid_decorrect(k);
         if(varmodel->weights[k]==inf) {
            par->val[k]=0;
         }
         else {
            collect_lhs_rhs(lhs, rhs, k);
            lhs += varmodel->weights[k];
            par->val[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         }
         resid_correct(k);
      }
      // need some thinking how to store sample info in file; likely every "saved" cycle.
      // main runs prepForOutput, and on a parVector main runs collectStats at the save intervals,
      // or it needs a new mechanism to switch on saving samples from output (which can be generic feature).
   }

   void sampleHpars() {
      varmodel->sample();
   }

   void restart() {
      varmodel->restart();
   }

   indepVarStr* varmodel;

};

// Here can start working on defining different variance structures for modelRreg
// For the moment, Rreg is only designed to accept variance models in the "indepVarStr" class, but this
// can include LASSO, BVS, DIAG / weighted, log-linear ...
class modelRregIden : public modelRreg {
public:
   modelRregIden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};

class modelRregDiag : public modelRreg {
public:
   modelRregDiag(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new diagVarStr(pmdescr, this->par);
   }
};

class modelRregGRL : public modelRreg {

   public:
   modelRregGRL(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod), beta_grid() {
      varmodel = new gridLVarStr(pmdescr, this->par);
      // make the grid, for the moment size and step are hard-coded in class member variabels below
      double sum_prob=0.0l;
      for(int i=0, double x=-grid_min_max; i<grid_size; i++, x+= grid_step) {
          grid_x[i] = x;
          grid_y[i] = exp(-abs(grid_x[i]));
          sum_prob += grid_y[i];
      }
      // scale grid-probs and store as log(prob's)
      for(int i=0; i<grid_size; i++) {
         grid_y[i] /= sum_prob;
         grid_y[i] = log(grid_y[i]);
      }
      beta_grid.initWith(M->ncol, 4); // 4 is middle value
 }

// modelRregGRL cannot use parent sample() and needs to re-define it
// some notes:
// - beta_grid is integer where each beta is in the grid, now from 0-8 and initialized at middle value 4
// - the 'x' axis values are grid_x[0-8], this now goes from -4 to 4 with step 1.0 (it is doubles because steps
//   are not neccessarily whole integers)
// - the 'y' axis values are grid_y[0-8], stored as logs
// - the actual beta-value includes the scale and is then beta[k] = scale*grid_x[beta_grid[k]]
// - the algorithm also uses difference of beta's (current - proposal), this can be directly computed from the
//   grid step-size and whether the proposal moves 'down' or 'up' in the grid (with move down current minus proposal
//   is positive step-size, with move up current minus proposal is negative step-size)
   void sample() {
      int curr_grid, prop_grid;
      double logtwo = log(2.0l);
      double loghalf = log(0.5l);
      double beta_diff;
      double lhs, rhs;
      double beta_scale = sqrt(varmodel->par[0]);
      double MHratio;
      for(size_t k=0; k < M->ncol; k++) {
         curr_grid = beta_grid[k];
         if(curr_grid == 0) {  // at left extreme, move up
            prop_grid = 1;
            beta_diff = -beta_scale*grid_step;
         }
         else if (curr_grid == 8) { // at right extreme, move down 
            prop_grid = 7;
            beta_diff = beta_scale*grid_step;
         }
         else {  // in between toss a coin how to move
            if(R::runif(0,1) < 0.5) {
               prop_grid = curr_grid-1;
               beta_diff = beta_scale*grid_step;
            }
            else {
               prop_grid = curr_grid + 1;
               beta_diff = -beta_scale*grid_step;
            }
         }
         collect_lhs_rhs(lhs, rhs, k);
         MHratio = -beta_diff*rhs - 0.5*beta_diff*beta_diff*lhs + 
            grid_prob[prop_grid] - grid_prob[curr_grid];
         if(curr_grid == 0 || curr_grid == 8) MHratio += loghalf;
         else if (prop_grid == 0 || prop_grid == 8) MHratio += logtwo;
         if(MHratio > 0 || log(R::runif(0,1)) < MHratio ) { // accept
            beta_grid[k] = prop_grid;
            resid_betaUpdate(beta_diff, k);
         }
      }
      // need to work on collecting fit that is needed for variance update
   }

   // Other methods from the parent class are:
   // - sampleHpars() - this calls varmodel->sample() and needs to have fit from this class - solve in constructor??
   // - restart() - this calls varmodel->restart() and is maybe not needed here? - insert an empty one in the variance class?

   int grid_size=9;
   double grid_min_max=4.0;
   double grid_step=1.0;
   double grid_x[9];
   double grid_y[9];
   simpleIntVector beta_grid;

}

class modelRregMixt : public modelRreg {
public:
   modelRregMixt(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new mixtVarStr(pmdescr, this->par);
   }
};

#endif /* modelRreg */
