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

#include "bootstrappers.h"


#include "gmwm_logic.h"

#include "guess_values.h"


// Include support functions
#include "transform_data.h"
#include "inline_functions.h"

// Functions that allow for manipulation
#include "armadillo_manipulations.h"

// Functions converted from R to Armadillo.
#include "rtoarmadillo.h"

// Functions that convert process to theoretical wv.
#include "process_to_wv.h"

// Functions that generate the process.
#include "gen_process.h"

// We use QMF, Haar, Select filter, etc.
#include "wv_filters.h"

// Wave Variance
#include "wave_variance.h"

// DWT
#include "dwt.h"

// Objective Functions
#include "objective_functions.h"

// rsample
#include "sampler.h"

// Covariance matrix
#include "covariance_matrix.h"

//' @title Bootstrap for Matrix V
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters, obtain the wavelet variance, and then V.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::mat cov_bootstrapper(const arma::vec&  theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                           unsigned int N, bool robust, double eff,
                           unsigned int H, bool diagonal_matrix){
  unsigned int nb_level = floor(log2(N));
  
  arma::mat res(nb_level, H);
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // MODWT transform
    arma::field<arma::vec> signal_modwt_bw = modwt_cpp(x, "haar", nb_level, "periodic", true);
    
    // Obtain WV
    arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
    // Store WV into matrix
    res.col(i) = wv_x;
  }
  
  // Do we need a diagnoal covariance matrix? 
  if(diagonal_matrix){
    return arma::diagmat(arma::cov(res.t()));
  }
  
  // Return the full covariance matrix
  return arma::cov(res.t());
}



//' @title Bootstrap for Optimism
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters, obtain the wavelet variance, and then V.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::mat optimism_bootstrapper(const arma::vec&  theta,
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                const arma::vec& scales, std::string model_type, 
                                unsigned int N, bool robust, double eff, double alpha,
                                unsigned int H){
  unsigned int nb_level = floor(log2(N));
  
  arma::mat theo(nb_level, H);
  
  arma::mat all_wv_empir(nb_level, H);
  
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // Obtain WV and confidence intervals
    arma::mat wvar = modwt_wvar_cpp(x, nb_level, robust, eff, alpha, "eta3", "haar");
    
    // Obtain the Omega matrix (CI HI, CI LO)
    arma::mat omega = arma::inv(fast_cov_cpp(wvar.col(2), wvar.col(1)));
    
    // Obtain the GMWM estimator's estimates. (WV_EMPIR)
    arma::vec est = gmwm_engine(theta, desc, objdesc, model_type, 
                                wvar.col(0), omega, scales, false);
    
    
    // Decomposition of the WV.
    theo.col(i) = theoretical_wv(est, desc, objdesc, scales);
    
    all_wv_empir.col(i) = wvar.col(0);
  }
  
  // Optimism Matrix bootstrap result 
  return cov(all_wv_empir.t(),theo.t());
}


//' @title Bootstrap for Optimism and GoF
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters, obtain the wavelet variance, and then V.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @param scales A \code{vec} containing the scales of the process.
//' @param model_type A \code{string} containing the model type either: SSM or IMU
//' @param N A \code{int} indicating how long the integer is. 
//' @param robust A \code{bool} indicating robust (T) or classical (F).
//' @param eff A \code{double} that handles efficiency.
//' @param H A \code{int} that indicates how many bootstraps should be obtained.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::field<arma::mat> opt_n_gof_bootstrapper(const arma::vec&  theta,
                                              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                              const arma::vec& scales, std::string model_type, 
                                              unsigned int N, bool robust, double eff, double alpha,
                                              unsigned int H){
  unsigned int nb_level = floor(log2(N));
  
  unsigned int p = theta.n_elem;
  
  arma::mat theo(nb_level, H);
  
  arma::mat all_wv_empir(nb_level, H);
  
  arma::vec obj_values(H);
  
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // Obtain WV and confidence intervals
    arma::mat wvar = modwt_wvar_cpp(x, nb_level, robust, eff, alpha, "eta3", "haar");
    
    // Obtain the Omega matrix (CI HI, CI LO)
    arma::mat omega = arma::inv(fast_cov_cpp(wvar.col(2), wvar.col(1)));
    
    // WV Empirical
    arma::vec wv_empir = wvar.col(0);
    
    // Take the mean of the first difference
    double expect_diff = mean_diff(x);
    
    // Min-Max / N
    double ranged = dr_slope(x);
    
    arma::vec theta_star = guess_initial(desc, objdesc, model_type, p, expect_diff, N, wvar, scales, ranged, 10000);
    
    // Obtain the GMWM estimator's estimates. (WV_EMPIR)
    arma::vec est = gmwm_engine(theta, desc, objdesc, model_type, 
                                wv_empir, omega, scales, false);
    
    arma::vec est_starting = gmwm_engine(theta_star, desc, objdesc, model_type, 
                                         wv_empir, omega, scales, true);
    
    // Obtain the objective value function
    obj_values(i) = getObjFun(est_starting, desc, objdesc, model_type, omega, wv_empir, scales); 
    
    // Decomposition of the WV.
    theo.col(i) = theoretical_wv(est, desc, objdesc, scales);
    
    all_wv_empir.col(i) = wvar.col(0);
  }
  
  // Out
  arma::field<arma::mat> out(2);
  
  // Optimism bootstrap result 
  out(0) = cov(all_wv_empir.t(),theo.t());
  
  // Obj Fun
  out(1) = obj_values; // Use N-1 and take by row
  
  // Optimism Matrix bootstrap result 
  return out;
}







//' @title Bootstrap for Standard Deviations of Theta Estimates
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::vec gmwm_sd_bootstrapper(const arma::vec&  theta,
                               const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                               const arma::vec& scales, std::string model_type,
                               unsigned int N, bool robust, double eff, double alpha,
                               unsigned int H){
  unsigned int nb_level = floor(log2(N));
  
  unsigned int p = theta.n_elem;
  
  arma::mat mest(H, p);
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // Obtain WV and confidence intervals
    arma::mat wvar = modwt_wvar_cpp(x, nb_level, robust, eff, alpha, "eta3", "haar");
    
    // Obtain the Omega matrix (CI HI, CI LO)
    arma::mat omega = arma::inv(fast_cov_cpp(wvar.col(2), wvar.col(1)));
    
    // Obtain the GMWM estimator's estimates. (WV_EMPIR)
    mest.col(i) = gmwm_engine(theta, desc, objdesc, model_type, 
             wvar.col(0), omega, scales, false);
    
  }
  
  // Return the sd of bootstrapped estimates
  return arma::stddev(mest.t());
}



//' @title Generate the Confidence Interval for GOF Bootstrapped
//' @description yaya
//' @param obj A \code{double} containing the objective value of the solution.
//' @param obj_boot A \code{vec} containing the objective values obtained while bootstrapping under Theta_hat.
//' @param B An \code{int} indicating how many times the bootstrapper should be run.
//' @param alpha A \code{double} indicating the amount of confidence for CI
//' @return A \code{vec} that has the alpha/2.0 quantile and then the 1-alpha/2.0 quantile. 
//' @keywords internal
// [[Rcpp::export]]
arma::vec boot_pval_gof(double obj, const arma::vec& obj_boot, unsigned int B = 1000, double alpha = 0.05){
  // Total number of objects in time series
  unsigned int N = obj_boot.n_elem;
  
  // Sequence of integers  
  arma::vec index = seq_cpp(0,N-1);
  
  arma::vec empty = arma::zeros<arma::vec>(0);
  
  // Bootstrap vector
  arma::vec boot(B);
  for (unsigned int b = 0; b < B; b++){
    arma::vec obj_boot_star = obj_boot.elem( arma::conv_to<arma::uvec>::from(rsample(index,N, true, empty)) );
    
    boot(b) = sum(obj < obj_boot_star)/double(N);
  } 
  
  arma::vec probs(2);
  probs(0) = alpha/2.0;
  probs(1) = 1.0-alpha/2;
  
  arma::vec quant = quantile_cpp(boot, probs);
  
  return quant;
}



//' @title Bootstrap for Estimating Both Theta and Theta SD
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters, obtain the wavelet variance, and then V.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::field<arma::mat> gmwm_param_bootstrapper(const arma::vec&  theta,
                                               const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                               const arma::vec& scales, std::string model_type,
                                               unsigned int N, bool robust, double eff, double alpha,
                                               unsigned int H){
  unsigned int nb_level = floor(log2(N));
  
  unsigned int p = theta.n_elem;
  
  arma::mat mest(p, H);
  
  arma::mat theo(nb_level, H);
  
  arma::mat wv_empir(nb_level, H);
  
  arma::vec obj_values(H);
  
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // Obtain WV and confidence intervals
    arma::mat wvar = modwt_wvar_cpp(x, nb_level, robust, eff, alpha, "eta3", "haar");
    
    // Obtain the Omega matrix (CI HI, CI LO)
    arma::mat omega = arma::inv(fast_cov_cpp(wvar.col(2), wvar.col(1)));
    
    // Obtain the GMWM estimator's estimates. (WV_EMPIR)
    mest.col(i) = gmwm_engine(theta, desc, objdesc, model_type, 
             wvar.col(0), omega, scales, false);
  }
  
  // Return the sd of bootstrapped estimates
  arma::field<arma::mat> out(2);
  
  // Theta estimate
  out(0) = mean(mest,1); // by row
  
  // Theta sd
  out(1) = stddev(mest,0,1); // Use N-1 and take by row
  
  return out;
}



//' @title Bootstrap for Everything!
//' @description Using the bootstrap approach, we simulate a model based on user supplied parameters, obtain the wavelet variance, and then V.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' Expand in detail...  
//' @author JJB
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::field<arma::mat> all_bootstrapper(const arma::vec&  theta,
                                        const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                        const arma::vec& scales, std::string model_type, 
                                        unsigned int N, bool robust, double eff, double alpha,
                                        unsigned int H){
  unsigned int nb_level = floor(log2(N));
  
  unsigned int p = theta.n_elem;
  
  arma::mat mest(p, H);
  
  arma::mat theo(nb_level, H);
  
  arma::mat all_wv_empir(nb_level, H);
  
  arma::vec obj_values(H);
  
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_model(N, theta, desc, objdesc);
    
    // Obtain WV and confidence intervals
    arma::mat wvar = modwt_wvar_cpp(x, nb_level, robust, eff, alpha, "eta3", "haar");
    
    // Obtain the Omega matrix (CI HI, CI LO)
    arma::mat omega = arma::inv(fast_cov_cpp(wvar.col(2), wvar.col(1)));
    
    arma::vec wv_empir = wvar.col(0);
    
    // Take the mean of the first difference
    double expect_diff = mean_diff(x);
    
    double ranged = dr_slope(x);
    
    // updated for WV obj pull
    arma::vec theta_star = guess_initial(desc, objdesc, model_type, p, expect_diff, N, wvar, scales, ranged, 10000);
    
    
    // Obtain the GMWM estimator's estimates. (WV_EMPIR)
    arma::vec est = gmwm_engine(theta, desc, objdesc, model_type, 
                                wv_empir, omega, scales, false);
    
    
    arma::vec est_starting = gmwm_engine(theta_star, desc, objdesc, model_type, 
                                         wv_empir, omega, scales, true);
    
    // Obtain the objective value function
    obj_values(i) = getObjFun(est_starting, desc, objdesc, model_type, omega, wv_empir, scales); 
    
    // Decomposition of the WV.
    theo.col(i) = theoretical_wv(est, desc, objdesc, scales);
    all_wv_empir.col(i) = wv_empir;
    // Store theta estimate
    mest.col(i) = est;
    
  }
  
  // Return the sd of bootstrapped estimates
  arma::field<arma::mat> out(5);
  
  all_wv_empir = all_wv_empir.t();
  
  // Optimism bootstrap result 
  out(0) = cov(all_wv_empir,theo.t());
  
  // V matrix
  out(1) = cov(all_wv_empir); // Take by row
  
  // Theta estimate
  out(2) = mean(mest,1); // by row
  
  // Theta sd
  out(3) = stddev(mest,0,1); // Use N-1 and take by row
  
  // Obj Fun
  out(4) = obj_values; // Use N-1 and take by row
  
  return out;
}