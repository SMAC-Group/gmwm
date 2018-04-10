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
#include <string>
// #include <omp.h>

#include "gmwm_logic.h"

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

// Covariance matrix
#include "covariance_matrix.h"

// Objective Functions
#include "objective_functions.h"

// Guess starting values
#include "guess_values.h"

// Count Models
#include "ts_checks.h"

// arima
#include "arima_gmwm.h"

//
#include "bootstrappers.h"


//' @title Optim loses NaN
//' @description This function takes numbers that are very small and sets them to the minimal tolerance for C++.
//' Doing this prevents NaN from entering the optim routine.
//' @param theta A \code{vec} that contains estimated GMWM theta values (untransformed).
//' @return A \code{vec} that contains safe theta values.
//' @keywords internal
// [[Rcpp::export]]
arma::vec code_zero(arma::vec theta){
  
  theta.elem(find(abs(theta) <= DBL_EPSILON)).fill(DBL_EPSILON);
  return theta;
}

//' @title Engine for obtaining the GMWM Estimator
//' @description This function uses the Generalized Method of Wavelet Moments (GMWM) to estimate the parameters of a time series model.
//' @param theta A \code{vec} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} containing a list of parameters (e.g. AR(1) = c(1,1), ARMA(p,q) = c(p,q,1))
//' @param model_type A \code{string} that represents the model transformation
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance
//' @param omega A \code{mat} that represents the covariance matrix.
//' @param scales A \code{vec} that contains the scales or taus (2^(1:J))
//' @param starting A \code{bool} that indicates whether we guessed starting (T) or the user supplied estimates (F).
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' If type = "imu" or "ssm", then parameter vector should indicate the characters of the models that compose the latent or state-space model.
//' The model options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"ARMA"}{an autoregressive moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' If model_type = "imu" or type = "ssm" then
//' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @backref src/gmwm_logic.cpp
//' @backref src/gmwm_logic.h
// [[Rcpp::export]]
arma::vec gmwm_engine(const arma::vec& theta,
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                      std::string model_type, 
                      arma::vec wv_empir,
                      arma::mat omega,
                      arma::vec scales,
                      bool starting){
  
  
  // Transform the Starting values
  arma::vec starting_theta = transform_values(theta, desc, objdesc, model_type);
                   
  // Apply Yannik's starting circle algorithm if our algorithm "guessed" the initial points
  if(starting){
    starting_theta = Rcpp_OptimStart(starting_theta, desc, objdesc, model_type, wv_empir, scales);
  }

  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
  
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, objdesc, model_type, omega, wv_empir, scales);
  
  return untransform_values(estim_GMWM, desc, objdesc, model_type);       

} 

//' @title Update Wrapper for the GMWM Estimator
//' @description This function uses information obtained previously (e.g. WV covariance matrix) to re-estimate a different model parameterization
//' @param theta A \code{vec} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} containing a list of parameters (e.g. AR(1) = c(1,1), ARMA(p,q) = c(p,q,1))
//' @param model_type A \code{string} that represents the model transformation
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance
//' @param omega A \code{mat} that represents the covariance matrix.
//' @param scales A \code{vec} that contains the scales or taus (2^(1:J))
//' @param starting A \code{bool} that indicates whether we guessed starting (T) or the user supplied estimates (F).
//' @return A \code{field<mat>} that contains the parameter estimates from GMWM estimator.
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @backref src/gmwm_logic.cpp
//' @backref src/gmwm_logic.h
// [[Rcpp::export]]
arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, unsigned int N, double expect_diff, double ranged, 
                                      const arma::mat& orgV, const arma::vec& scales, const arma::mat& wv,
                                      bool starting, 
                                      std::string compute_v, unsigned int K, unsigned int H,
                                      unsigned int G, 
                                      bool robust, double eff){
  
  // Number of parameters
  unsigned int np = theta.n_elem;
    
  // Guessed Values
  arma::vec guessed_theta = theta;
  
  // V matrix
  arma::mat V = orgV;
  
  // Diagonal Omega Matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  arma::vec wv_empir = wv.col(0);
  
  // Do we need to run a guessing algorithm?
  if(starting){

    theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wv, scales, ranged, G);
    
    guessed_theta = theta;
  }

  // Obtain the GMWM estimator estimates.
  theta = gmwm_engine(theta, desc, objdesc, model_type, 
                      wv_empir, omega, scales, starting);

  theta = code_zero(theta);
    
  // Bootstrap the V matrix
  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        // False here means we create the "full V" matrix
        V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false);
        omega = arma::inv(diagmat(V));
        
        // The theta update in this case MUST not use Yannick's starting algorithm. Hence, the false value.
        theta = gmwm_engine(theta, desc, objdesc, model_type, wv_empir, omega, scales, false);
        
        // Optim may return a very small value. In this case, instead of saying its zero (yielding a transform issue), make it EPSILON.
        theta = code_zero(theta);
    }
  }
  
  
  std::map<std::string, int> models = count_models(desc);
  
  // Order AR1s so largest phi is first!
  if(models["AR1"] > 1  || models["GM"] > 1){
    theta = order_AR1s(theta, desc, objdesc);
  }
  
  // Obtain the theoretical WV.
  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);
  
  // Obtain the objective value function
  arma::vec obj_value(1);
  obj_value(0) = getObjFun(theta, desc, objdesc,  model_type, omega, wv_empir, scales); 
  
  // Export calculations to R.
  arma::field<arma::mat> out(6);
  out(0) = theta;
  out(1) = guessed_theta;
  out(2) = V;
  out(3) = theo;
  out(4) = decomp_theo;
  out(5) = obj_value;
  
  return out;
                                        
}

//' @title Master Wrapper for the GMWM Estimator
//' @description This function generates WV, GMWM Estimator, and an initial test estimate.
//' @param data A \code{vec} containing the data.
//' @param theta A \code{vec} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param objdesc A \code{field<vec>} containing a list of parameters (e.g. AR(1) = c(1,1), ARMA(p,q) = c(p,q,1))
//' @param model_type A \code{string} that represents the model transformation
//' @param starting A \code{bool} that indicates whether the supplied values are guessed (T) or are user-based (F).
//' @param alpha A \code{double} that handles the alpha level of the confidence interval (1-alpha)*100
//' @param compute_v A \code{string} that describes what kind of covariance matrix should be computed.
//' @param K An \code{int} that controls how many times theta is updated.
//' @param H An \code{int} that controls how many bootstrap replications are done.
//' @param G An \code{int} that controls how many guesses at different parameters are made.
//' @param robust A \code{bool} that indicates whether the estimation should be robust or not.
//' @param eff A \code{double} that specifies the amount of efficiency required by the robust estimator.
//' @return A \code{field<mat>} that contains a list of ever-changing estimates...
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @export
//' @backref src/gmwm_logic.cpp
//' @backref src/gmwm_logic.h
// [[Rcpp::export]]
arma::field<arma::mat> gmwm_master_cpp(arma::vec& data, 
                                      arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, bool starting,
                                      double alpha, 
                                      std::string compute_v, unsigned int K, unsigned int H,
                                      unsigned int G, 
                                      bool robust, double eff){
  
  // Obtain counts of the different models we need to work with
  std::map<std::string, int> models = count_models(desc);
  
  // HACK METHOD (todo: formalize it)
  // Determine if we need to difference
  if(models["SARIMA"] > 0){
    
    // Note: s, i, si are 6,7,8 => 5,6,7
    for(unsigned int i = 0; i < desc.size(); i ++){
      if(objdesc(i).n_elem > 3){
        arma::vec sarima_desc = objdesc(i);
        // Do we need to difference? 
        if(sarima_desc(6) > 0 || sarima_desc(7) > 0){
          // Perform differencing in specific order...
          // First non-seasonal and then seasonal.
          
          if(sarima_desc(6) > 0){
            // Lag is always 1, number of differences is (i)
            data = diff_cpp(data, 1, sarima_desc(6));
          }
          
          if(sarima_desc(7) > 0){
            // Lag is always seasonality (s) and differences is (si).
            data = diff_cpp(data, sarima_desc(5), sarima_desc(7));
          }
          
          // Kill loop. We only handle the tsmodel object with the first difference.
          break;
        }
      }
    }
  }
  
  // ------ Variable Declarations
  
  // Length of the Time Series
  unsigned int N = data.n_elem;
  
  // Number of Scales (J)
  unsigned int nlevels = floor(log2(N));
  
  // Number of parameters
  unsigned int np = theta.n_elem;
  
  // Take the mean of the first difference
  double expect_diff = mean_diff(data);
  
  // Guessed values of Theta (user supplied or generated)
  arma::vec guessed_theta = theta;
  
  // MODWT decomp
  arma::field<arma::vec> modwt_decomp = modwt_cpp(data, "haar", nlevels, "periodic", true);
  
  // Obtain WV and confidence intervals
  arma::mat wvar = wvar_cpp(modwt_decomp, robust, eff, alpha, "eta3");
  
  // Extract
  arma::vec wv_empir = wvar.col(0);
  arma::vec ci_lo = wvar.col(1);
  arma::vec ci_hi = wvar.col(2);
  
  //-------------------------
  // Obtain Covariance Matrix
  //-------------------------
  
  arma::mat V;
  
  // compute_cov_cpp is the hard core function. It can only be improved by using parallelization.
  if(compute_v == "diag" || compute_v == "full"){
    arma::field<arma::mat> Vout = compute_cov_cpp(modwt_decomp, nlevels, compute_v, robust, eff);
    if(robust){
      V = Vout(1);
    }else{
      V = Vout(0);
    }
  }else{
     V = fast_cov_cpp(wvar.col(2), wvar.col(1));
  }

  // Obtain the Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Store the original V matrix (in case of bootstrapping) for use in the update function
  arma::mat orgV = V;
  
  // Calculate the values of the Scales 
  arma::vec scales = scales_cpp(nlevels);
  
  // Min-Max / N
  double ranged = dr_slope(data);
  
  // Guess starting values for the theta parameters
  if(starting){
    
    // Always run guessing algorithm
    theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wvar, scales, ranged, G);
    
    // If under ARMA case and only ARMA is in the model, 
    // then see how well these values are.
    if(desc.size() == 1 && (desc[0] == "SARIMA" || desc[0] == "ARMA11") && N <= 1000){
      
      // Use R's ARIMA function to estimate parameter space
      arma::vec theta2 = Rcpp_ARIMA(data, objdesc(0)); // Only 1 objdesc in the available.
      
      // Obtain the obj function under omega with these initial guesses
      // DO >>NOT<< USE Yannick's to optimize starting values!!!!
      double mle_css_obj = getObjFun(theta2, desc, objdesc,  model_type, omega, wv_empir, scales); 
      
      // Obtain the objective function under Yannick's starting algorithm
      double init_guess_obj = getObjFunStarting(theta, desc, objdesc,  model_type, wv_empir, scales);
      
      // What performs better? 
      if(mle_css_obj < init_guess_obj){
        // Disable starting value optimization if using MLE. 
        theta = theta2;
        starting = false;
      }

    }
    
    guessed_theta = theta;
  }
  
  // Obtain the GMWM estimator's estimates.
  theta = gmwm_engine(theta, desc, objdesc, model_type, 
                      wv_empir, omega, scales, starting);
  
  // Optim may return a very small value. In this case, instead of saying its zero (yielding a transform issue), make it EPSILON.
  theta = code_zero(theta);
  
  // Enable bootstrapping
  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        // Create the full V matrix
        V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false);
      
        // Update the omega matrix
        omega = arma::inv(diagmat(V));
        
        // The theta update in this case MUST not use Yannick's starting algorithm. Hence, the false value.
        theta = gmwm_engine(theta, desc, objdesc, model_type, wv_empir, omega, scales, false);
        
        // Optim may return a very small value. In this case, instead of saying its zero (yielding a transform issue), make it EPSILON.
        theta = code_zero(theta);
    }
  }

  
  if(desc[0] == "SARIMA" && desc.size() == 1){
    
    arma::vec temp = objdesc(0);
    unsigned int p = temp(0);
    if(p != 0 && invert_check(arma::join_cols(arma::ones<arma::vec>(1), -theta.rows(0, p - 1))) == false){
      Rcpp::Rcout << "WARNING: This ARMA model contains AR coefficients that are NON-STATIONARY!" << std::endl;
    }
  } 
  
  // Order AR1s / GM so largest phi is first!
   if(models["AR1"] > 1 || models["GM"] > 1){
     theta = order_AR1s(theta, desc, objdesc);
   }

  // Obtain the objective value function
  arma::vec obj_value(1);
  obj_value(0) = getObjFun(theta, desc, objdesc,  model_type, omega, wv_empir, scales); 
  
  arma::vec dr_s(1);
  dr_s(0) = ranged;
  
  // Decomposition of the WV.
  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);

  // Export information back
  arma::field<arma::mat> out(13);
  out(0) = theta;
  out(1) = guessed_theta;
  out(2) = wv_empir;
  out(3) = ci_lo;
  out(4) = ci_hi;
  out(5) = V;
  out(6) = orgV;
  out(7) = expect_diff;
  out(8) = theo;
  out(9) = decomp_theo;
  out(10) = obj_value;
  out(11) = omega;
  out(12) = dr_s;
  return out;
}




/*

//' @title Simulate GMWM
//' 
//' @examples
//' x=rnorm(100)
//' wavelet_variance_cpp(x, "haar", "diag")
// [[Rcpp::export]]
arma::field<arma::mat> simGMWM(const arma::vec& theta, const arma::mat& omega,
                               const std::vector<std::string>& desc, const arma::vec& wv_empir,
                               const arma::vec& tau, unsigned int N, unsigned int B = 500, bool var_or_mu = false){
  
  // Number of parameters
  unsigned int num_param = theta.n_elem;
  
  // Initialisation of results structures
  arma::mat GMWM(B,num_param);
  arma::mat GMWM_plus(B,num_param);
  
  // Starting values
	arma::vec starting_theta = set_result_values(theta, desc);

  // Start bootstrap
  for(unsigned int b=0; b<B; b++){  	
  	// Simulate  	
  	arma::vec x = gen_model(N, theta, desc);
    
  	// ------------------------------------
  	// Compute standard GMWM
  	// ------------------------------------
  	// Compute WV
  	arma::field<arma::mat> wv_x = wavelet_variance_cpp(x, "haar", "diag");
  	
  	// Omega matrix
  	arma::mat omega = arma::inv(diagmat(wv_x(4)));
  	
  	// Empirical WV
  	arma::vec wv_empir = wv_x(0);
      	
  	// Find GMWM estimator
  	arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, omega, tau, wv_empir, N);
  	
  	// Save results
  	GMWM.row(b) = set_result_values(estim_GMWM, desc);
  	
  	// ------------------------------------
  	// Compute augmented GMWM
  	// ------------------------------------
  	// Compute Omega
  	arma::mat V = cov_bootstrapper(GMWM.row(b), desc, max(tau), N, B, var_or_mu);
  	arma::mat omega_v = arma::inv(diagmat(V));
  	
  	// Empirical WV + variance
    arma::vec temp(1);
    if(var_or_mu){
      temp(0) = arma::var(x);
    }
    else{
      temp(0) = arma::mean(x);
    }
    
  	arma::vec wv_empir_v = join_cols(temp,wv_empir);
  	
  	// Find GMWM+ estimator
  	arma::vec estim_GMWM_plus = Rcpp_Optim(theta, desc, omega_v, tau, wv_empir_v, N);
  	
  	// Save results
  	GMWM_plus.row(b) = set_result_values(estim_GMWM_plus, desc);
  }
  
  arma::field<arma::mat> out(2);
  out(0) = GMWM;
  out(1) = GMWM_plus;
  
  return out;
} 


//
// Mondal and Percival estimator
//

//What/Where is get.slepians?
arma::vec percival(arma::vec x){
  arma::vec xsq = arma::square(x);
  double Tn = log(median(xsq));
  arma::mat beta = get.slepians(npoints=x.n_elem,nwin=5)/x.n_elem;
  arma::rowvec colsums = arma::sum(beta); // 1 x n
  arma::vec J = arma::trans(beta)*arma::sign(log(xsq)-Tn);
  arma::vec mu = (J * colsums)/(colsums*arma::trans(colsums));
  arma::mat Ahat = arma::mean(arma::square(J-mu*colsums));
  
  double temp = R::qnorm(3.0/4.0, 0.0, 1.0, 1,0); //0.6744898
  // dnorm(qnorm(3/4)) = 0.3177766
  
  arma::vec muhat = Tn-2*log(temp)-(Ahat/square(-2*R::dnorm(temp, 0.0, 1.0, 0)*temp))/(2*x.n_elem);
  return exp(muhat);
}*/