/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
#include "robust_components.h"
#include "inline_functions.h"
using namespace Rcpp;

/* -------------------------------- Start ROBUST FUNCTIONS -------------------------- */

// @title Objective Function for Tuning Constant
// @description Objective function that finds tuning constant
// @usage objFun_find_biwc(crob, eff)
// @param crob A \code{double} that indicates the value to optimize. Normally between [0,5].
// @param eff A \code{double} that represents the efficiency.
// @return A \code{double} that contains the biwc.
// @details This function is not available from within R.
// @name objFun_find_biwc
// @docType methods
// @rdname objFun_find_biwc-methods
double objFun_find_biwc(double crob, double eff){
  
  // q, mean, sd, lower.tail, log.p
  double norm_prob = R::pnorm(crob,0.0,1.0,1,0);
  double norm_density = R::dnorm(crob,0.0,1.0,0);
  
  // avoid using pow as much as possible
  double crob_sq = square(crob);
  double crob_quad = square(crob_sq);
  double crob_eight = square(crob_quad);
  
  double mu20c=1309458150*norm_prob-2*crob*(654729075+218243025*crob_sq+43648605*crob_quad+6235515*(crob_sq*crob_quad)+692835*crob_eight+62985*(crob_eight*crob_sq)+4845*(crob_quad*crob_eight)+323*(crob_sq*crob_quad*crob_eight)+19*(crob_eight*crob_eight)+(crob_eight*crob_eight*crob_sq))*norm_density-654729075;
  double mu18c=68918850*norm_prob-2*crob*(34459425+11486475*crob_sq+2297295*crob_quad+328185*(crob_sq*crob_quad)+36465*crob_eight+3315*(crob_eight*crob_sq)+255*(crob_quad*crob_eight)+17*(crob_sq*crob_quad*crob_eight)+(crob_eight*crob_eight))*norm_density-34459425;
  double mu16c=4054050*norm_prob-2*crob*(2027025+675675*crob_sq+135135*crob_quad+19305*(crob_sq*crob_quad)+2145*crob_eight+195*(crob_eight*crob_sq)+15*(crob_quad*crob_eight)+(crob_sq*crob_quad*crob_eight))*norm_density-2027025;
  double mu14c=270270*norm_prob-2*crob*(135135+45045*crob_sq+9009*crob_quad+1287*(crob_sq*crob_quad)+143*crob_eight+13*(crob_eight*crob_sq)+(crob_quad*crob_eight))*norm_density-135135;
  double mu12c=20790*norm_prob-2*crob*(10395+3465*crob_sq+693*crob_quad+99*(crob_sq*crob_quad)+11*crob_eight+(crob_eight*crob_sq))*norm_density-10395;
  double mu10c=1890*norm_prob-2*crob*(945+315*crob_sq+63*crob_quad+9*(crob_sq*crob_quad)+crob_eight)*norm_density-945;
  double mu8c=210*norm_prob-2*crob*(105+35*crob_sq+7*crob_quad+(crob_sq*crob_quad))*norm_density-105;
  double mu6c=30*norm_prob-2*crob*(15+5*crob_sq+crob_quad)*norm_density-15;
  double mu4c=6*norm_prob-2*crob*(3+crob_sq)*norm_density-3;
  double mu2c=2*norm_prob-2*crob*norm_density-1;
  double ac=(1/crob_eight)*mu10c-(4/(crob_sq*crob_quad))*mu8c+(6/crob_quad)*mu6c-(4/crob_sq)*mu4c+mu2c;
  double Q=(28/crob_quad)*mu8c-(8/crob_sq)*mu6c+(1/(crob_eight*crob_eight))*mu20c+(8/(crob_sq*crob_quad*crob_eight))*mu18c+(28/(crob_quad*crob_eight))*mu16c-(56/(crob_eight*crob_sq))*mu14c+(70/crob_eight)*mu12c-(56/(crob_sq*crob_quad))*mu10c+mu4c-ac*ac;
  double M=(1/crob_eight)*mu12c-(4/(crob_sq*crob_quad))*mu10c+(6/crob_quad)*mu8c-(4/crob_sq)*mu6c+mu4c-ac;
  return square((0.5*M*M)/Q - eff);
}

// @title Obtain Tuning Constant crob.bw
// @description Objective function that finds tuning constant
// @usage find_biwc(eff)
// @param eff A \code{double} that represents the desired level of efficiency as compared to the classic MLE.
// @return A \code{double} that contains the optimized biwc.
// @examples
// find_biwc(0.6)
double find_biwc(double eff = 0.6){
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optimize = stats["optimize"];    

  Rcpp::List crob_t = optimize(_["f"]  = Rcpp::InternalFunction(&objFun_find_biwc),
                          _["lower"] = 0,
                          _["upper"] = 100,
                          _["eff"] = eff);
                          
  double out = as<double>(crob_t[0]);

  return out;
}


// Objective funcion for robust WV
double objFun_sig_rob_bw(double sig2_bw, arma::vec x, double a_of_c, double crob_bw){
  arma::vec r = x/sqrt(sig2_bw);
  arma::vec rsq = arma::square(r);
  arma::uvec ivec = (abs(r) > crob_bw);
  arma::vec w = ((1 - ivec) % arma::square(1 - rsq/square(crob_bw)) );
  return square(arma::mean(arma::square(r)%arma::square(w)) - a_of_c);
}

// Robust estimator. Inputs are wavelet coefficients (y) and desired level of efficiency. 
double sig_rob_bw(arma::vec y, double eff = 0.6){
  double crob_bw = find_biwc(eff);

  arma::vec x = y/arma::stddev(y);  
  
  // q, mean, sd, lower.tail, log.p
  double norm_prob = R::pnorm(crob_bw,0.0,1.0,1,0);
  double norm_density = R::dnorm(crob_bw,0.0,1.0,0);
  
  // avoid using pow as much as possible
  double crob_sq = square(crob_bw);
  double crob_quad = square(crob_sq);
  double crob_eight = square(crob_quad);
  
  double a_of_c = (1/crob_eight)*(1890*norm_prob-2*crob_bw*(945+315*crob_sq+63*crob_quad+9*(crob_sq*crob_quad)+crob_eight)*norm_density-945)
                  -(4/(crob_sq*crob_quad))*(210*norm_prob-2*crob_bw*(105+35*crob_sq+7*crob_quad+(crob_sq*crob_quad))*norm_density-105)
                  +(6/crob_quad)*(30*norm_prob-2*crob_bw*(15+5*crob_sq+crob_quad)*norm_density-15)
                  -(4/crob_sq)*(6*norm_prob-2*crob_bw*(3+crob_sq)*norm_density-3)
                  +2*norm_prob-2*crob_bw*norm_density-1;

  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optimize = stats["optimize"];    

  Rcpp::List opt_val = optimize(_["f"]  = Rcpp::InternalFunction(&objFun_sig_rob_bw),
                          _["lower"] = 0,
                          _["upper"] = 2,
                          _["x"] = x,
                          _["a_of_c"] = a_of_c,
                          _["crob_bw"] = crob_bw);
          
  double sig2_hat_rob_bw = as<double>(opt_val[0])*var(y);
  return sig2_hat_rob_bw;
}

/* -------------------------------- End ROBUST FUNCTIONS -------------------------- */