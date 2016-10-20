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

//#include "summary_inf_mod.h"

// Inference goodies
#include "inference.h"

// Model Selection criterion
#include "analytical_matrix_derivatives.h"

// Bootstraps!
#include "bootstrappers.h"


//' @title Routing function for summary info
//' @description Gets all the data for the summary.gmwm function.
//' @param theta A \code{vec} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} containing a list of parameters (e.g. AR(1) = c(1,1), ARMA(p,q) = c(p,q,1))
//' @param model_type A \code{string} that represents the model transformation
//' @param wv_empir A \code{vec} that 
//' @param theo A \code{vec} that
//' @param scales A \code{vec} that
//' @param V A \code{mat} that contains the V matrix used to obtain the GMWM.
//' @param omega A \code{mat} that 
//' @param obj_value A \code{double} that contains the objective function value at the optimized solution.
//' @param N A \code{int} that indicates how long the time series is.
//' @param alpha A \code{double} that handles the alpha level of the confidence interval (1-alpha)*100
//' @param robust A \code{bool} that indicates whether the estimation should be robust or not.
//' @param eff A \code{double} that specifies the amount of efficiency required by the robust estimator.
//' @param inference A \code{bool} that indicates whether inference (e.g. GoF) should be run.
//' @param fullV A \code{bool} that indicates whether the matrix has been fully bootstrapped.
//' @param bs_gof A \code{bool} indicating whether the GoF should be bootstrapped or done asymptotically.
//' @param bs_gof_p_ci A \code{bool} indicating whether a bootstrapped p-value should be generated during the bootstrapped GoF
//' @param bs_ci A \code{bool} that indicates whether a bootstrapped CI should be obtained or to use analytical derivatives.
//' @param B A \code{int} that indicates how many iterations should take place.
//' @return A \code{field<mat>} that contains bootstrapped / asymptotic GoF results as well as CIs.
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::mat> get_summary(arma::vec theta,
                                   const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                   std::string model_type, 
                                   const arma::vec& wv_empir, const arma::vec& theo, const arma::vec& scales,
                                   arma::mat V, const arma::mat& omega, double obj_value,
                                   unsigned int N, double alpha,
                                   bool robust, double eff, 
                                   bool inference, bool fullV,
                                   bool bs_gof,  bool bs_gof_p_ci, bool bs_theta_est, bool bs_ci, 
                                   unsigned int B){
  
  // Inference test result storage
  arma::vec gof;
  
  arma::mat ci;
  

  // Bootstrap results
  arma::vec sd;

  arma::vec bs_obj_values;
    

  // Determine the type of bootstrapper needed.
  
  if(!(bs_ci || bs_gof) & !fullV & inference){
     V = cov_bootstrapper(theta,
                       desc, objdesc,
                       N, robust, eff,
                       B, true);
  }else if(bs_ci || bs_gof){
    arma::field<arma::mat> bs = all_bootstrapper(theta,
                                                 desc, objdesc,
                                                 scales, model_type, 
                                                 N, robust, eff, alpha, B);
    
    if(!fullV){
      V = bs(1);
    }
    
    if(bs_theta_est){
      theta = bs(2);
    }
    
    if(bs_ci){
      sd = bs(3);
    }
    
    if(bs_gof){
      bs_obj_values = bs(4);
    }
  }

  // Generate inference information
  if(inference){
    // Obtain a confidence interval for the parameter estimates AND calculate chisq goodness of fit
    if(bs_ci){
      ci = format_ci(theta, sd, alpha);
    }else{
      arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);
      ci = theta_ci(theta, A, V, omega, alpha);
    }

    if(bs_gof){
      gof = bootstrap_gof_test(obj_value, bs_obj_values, alpha, bs_gof_p_ci); 
    }else{
      gof = gof_test(theta, desc, objdesc, model_type, scales, V, wv_empir);
    }

  }


  // Export information back
  arma::field<arma::mat> out(2);
  out(0) = ci;
  out(1) = gof;
  
  return out;
  
}