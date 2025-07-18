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
#include <cmath>
#include "modelMatrix.h"
#include "indepVarStr.h"
#include "dataMatrix.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "modelHelper.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelMatrix(modeldescr, rmod)   {
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

class modelRregGRL : public modelRreg {                                      // *** GRid Lasso ***

   public:
   modelRregGRL(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod), beta_grid() {
      varmodel = new gridLVarStr(pmdescr, this->par);
      beta_grid.initWith(M->ncol, grid.mid); // initialize grid-steps as the middle value
      ppi = new modelHelper(pmdescr, 0.0l, *(this->par), "ppi");
      // it is important to start the gridLASSO at a roughly right scale, it is taken here
      // as 0.10 * (raw response var) / Npredictors
      varmodel->par->val[0]= 0.1*rmod->stats.var/double(M->ncol);
   }

   ~modelRregGRL() {
      delete ppi;    // the varmodel is deleted in the parent
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
      double beta_scale = sqrt(varmodel->par->val[0]);
      double MHratio;
//      int count_accept=0;
      for(size_t k=0; k < M->ncol; k++) {
         curr_grid = beta_grid[k];
         if(curr_grid == 0) prop_grid = 1;                           // at left extreme, move up
         else if (curr_grid == grid.last) prop_grid = grid.last - 1; // at right extreme, move down 
         else {                                                      // in between toss a coin how to move
            if(R::runif(0,1) < 0.5) prop_grid = curr_grid-1;
            else prop_grid = curr_grid + 1;
         }
         beta_diff = beta_scale*(grid.x[curr_grid]-grid.x[prop_grid]);
         collect_lhs_rhs(lhs, rhs, k);
         MHratio = -beta_diff*rhs - 0.5*beta_diff*beta_diff*lhs + 
            grid.logp[prop_grid] - grid.logp[curr_grid];
         if(curr_grid == 0 || curr_grid == grid.last) MHratio += loghalf;
         else if (prop_grid == 0 || prop_grid == grid.last) MHratio += logtwo;
         if(MHratio > 0 || log(R::runif(0,1)) < MHratio ) { // accept
            beta_grid[k] = prop_grid;
            par->val[k] = beta_scale*grid.x[prop_grid];
            resid_fit_betaUpdate(beta_diff, k);
//            count_accept++;
           // also update ppi, it is 0 when at the mode, and 1 otherwise
           if(beta_grid[k]==grid.mid) ppi->par->val[k] = 0.0l;
           else ppi->par->val[k] = 1.0l;
         }
      }
   }

   // also sampleHpars needs to be re-defined, the version in the parent class calls the usual
   // varmodel->sample(), but here the varmodel->sampleScale() should be used.
   void sampleHpars() {
      double oldscale = sqrt(varmodel->par->val[0]);
      double lhs, rhs;
      getFitScaleStats(lhs, rhs);
      varmodel->sampleScale(lhs, rhs);
      resid_fit_scaleUpdate(oldscale,sqrt(varmodel->par->val[0]));
   }

   // Standard Bayesian LASSO-based grid
   /*
   struct {size_t n {7}; size_t mid {3}; size_t last {6};
            double x[7] {-5, -3.5, -2, 0, 2, 3.5, 5}; 
            double p[7] {0.005, 0.022, 0.101, 0.744, 0.101, 0.022, 0.005};
            double logp[7] {-5.296, -3.796, -2.296, -0.296, -2.296, -3.796, -5.296}; } grid;
   */

   // Epow(0.7) grid
   /*
   struct {size_t n {7}; size_t mid {3}; size_t last {6};
   double x[7] {-10, -6, -3, 0, 3, 6, 10}; 
   double p[7] {0.005, 0.023, 0.089, 0.767, 0.089, 0.023, 0.005};
   double logp[7] {-5.278, -3.771, -2.424, -0.266, -2.424, -3.771, -5.278}; } grid;
   */

   // Epow(0.5) grid
   struct {size_t n {7}; size_t mid {3}; size_t last {6};
            double x[7] {-20, -10, -5, 0, 5, 10, 20}; 
            double p[7] {0.009, 0.032, 0.081, 0.757, 0.081, 0.032, 0.009};
            double logp[7] {-4.751, -3.441, -2.515, -0.279, -2.515, -3.441, -4.751}; } grid;

//   int grid_size=9;
//   double grid_min_max=4.0;
//   double grid_step=1.0;
//   double grid_x[9];
//   double grid_y[9];

   simpleIntVector beta_grid;
   modelHelper* ppi;

};

class modelRregMixt : public modelRreg {
public:
   modelRregMixt(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new mixtVarStr(pmdescr, this->par);
   }
};

#endif /* modelRreg */
