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

#include "guess_values.h"

using namespace arma;
using namespace Rcpp;


//' @title Guided Initial Values for GMWM Estimator with Starting Technique
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param V A \code{matrix} that represents the covariance matrix.
//' @param wv_empir A \code{vector} that contains the empirical wavelet variance
//' @param N A \code{integer} that indicates the length of the signal being studied.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). If type = "IMU" or "SSM", then
//' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
//' options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"ARMA"}{an autoregressiveß moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' If type = "ARMA", the function takes condition least squares as starting values; if type = "IMU" or type = "SSM" then
//' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec gmwm_cpp(const arma::vec& theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type, 
                          const arma::mat& V, const arma::vec& wv_empir,
                          const arma::vec& tau){
                                 
  // Number of parameters
  //unsigned int num_param = theta.n_elem;
    
  // Starting values
  arma::vec starting_theta = transform_values(theta, desc, objdesc, model_type);
  
  // Optimize Starting values via Jannick's Method
  starting_theta = Rcpp_OptimStart(starting_theta, desc, objdesc, model_type, wv_empir, tau);
  
  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
      
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, objdesc, model_type, omega, wv_empir, tau);

  return trans(untransform_values(estim_GMWM, desc, objdesc, model_type));                          
}

//' @title User Specified Initial Values for GMWM Estimator
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model.
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param V A \code{matrix} that represents the covariance matrix.
//' @param wv_empir A \code{vector} that contains the empirical wavelet variance
//' @param N A \code{integer} that indicates the length of the signal being studied.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). If type = "IMU" or "SSM", then
//' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
//' options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"ARMA"}{an autoregressiveß moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' If type = "ARMA", the function takes condition least squares as starting values; if type = "IMU" or type = "SSM" then
//' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec adv_gmwm_cpp(const arma::vec& theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type, 
                          const arma::mat& V, const arma::vec& wv_empir,
                          const arma::vec& tau){

  // Starting values
  arma::vec starting_theta = transform_values(theta, desc, objdesc, model_type);
    
  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
      
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, objdesc, model_type, omega, wv_empir, tau);

  return trans(untransform_values(estim_GMWM, desc, objdesc, model_type));                          
}


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
arma::mat gmwm_bootstrapper(const arma::vec&  theta,
                            const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                            unsigned int N, bool robust, double eff,
                            unsigned int H){
  unsigned int nb_level = floor(log2(N));
    
  arma::mat res(H, nb_level);
	for(unsigned int i=0; i<H; i++){
		arma::vec x = gen_model(N, theta, desc, objdesc);

    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level, "periodic");
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter(), "modwt");
  
		arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
	  res.row(i) = arma::trans(wv_x);
	}
	return arma::diagmat(arma::cov(res));
}

// [[Rcpp::export]]
arma::vec gmwm_engine(const arma::vec& theta,
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                      std::string model_type,  bool robust,
                      arma::vec wv_empir,
                      arma::mat V,
                      arma::vec scales,
                      bool starting){
  /*unsigned int np = theta.n_elem;  
  if(robust){
    // Params np
    np = np + 1;
    
    scales = scales.rows(0,np-1);
    wv_empir = wv_empir.rows(0,np-1);
    // .submat( first_row, first_col, last_row, last_col )
    V = V.submat(0, 0, np-1, np-1 );
    
    np = np - 1;
  }*/
  
  // Apply Yannik's starting circle algorithm if we "guessed" the initial points or a user overrides it
  arma::rowvec estimate;
  if(starting){
    estimate = gmwm_cpp(theta, desc, objdesc, model_type, V, wv_empir, scales);
  }else{
    estimate = adv_gmwm_cpp(theta, desc, objdesc, model_type, V, wv_empir, scales);
  }

  return trans(estimate);  
} 

// [[Rcpp::export]]
arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, unsigned int N, double expect_diff, 
                                      arma::mat orgV, arma::vec scales, arma::vec wv_empir,
                                      bool starting, 
                                      std::string compute_v, unsigned int K, unsigned int H,
                                      unsigned int G, 
                                      bool robust, double eff){
  
  unsigned int np = theta.n_elem;
    
  arma::vec guessed_theta = theta;
  
  arma::mat V = orgV;
    
  if(starting){

    theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wv_empir, scales, G);
    
    guessed_theta = theta;
  }

  
  theta = gmwm_engine(theta, desc, objdesc, model_type, robust,
                      wv_empir, V, scales, starting);
  //

  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        V = gmwm_bootstrapper(theta, desc, objdesc, N, robust, eff, H);
        theta = gmwm_engine(theta, desc, objdesc, model_type, robust,
                      wv_empir, V, scales, starting);
    }
  }

  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);

  arma::field<arma::mat> out(5);
  out(0) = theta;
  out(1) = guessed_theta;
  out(2) = V;
  out(3) = theo;
  out(4) = decomp_theo;
  
  return out;
                                        
}

// [[Rcpp::export]]
arma::field<arma::mat> gmwm_master_cpp(const arma::vec& data, 
                                      arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, bool starting,
                                      double p, 
                                      std::string compute_v, unsigned int K, unsigned int H,
                                      unsigned int G, 
                                      bool robust, double eff){
  unsigned int N = data.n_elem;
  unsigned int nlevels = floor(log2(N));
  
  unsigned int np = theta.n_elem;
  
  double expect_diff = mean_diff(data);
  
  arma::vec guessed_theta = theta;
  
  arma::field<arma::vec> modwt_decomp = modwt_cpp(data, "haar", nlevels, "periodic");
  
  arma::mat wvar = wvar_cpp(modwt_decomp, robust, eff, p, "eta3", "haar");
  
  arma::vec wv_empir = wvar.col(0);
  arma::vec ci_lo = wvar.col(1);
  arma::vec ci_hi = wvar.col(2);
  
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
  
  arma::mat orgV = V;
  
  arma::vec scales = scales_cpp(nlevels);
  
  if(starting){

    theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wv_empir, scales, G);
    
    guessed_theta = theta;
  }

  
  theta = gmwm_engine(theta, desc, objdesc, model_type, robust,
                      wv_empir, V, scales, starting);
  //

  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        V = gmwm_bootstrapper(theta, desc, objdesc, N, robust, eff, H);
        theta = gmwm_engine(theta, desc, objdesc, model_type, robust,
                      wv_empir, V, scales, starting);
    }
  }

  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);

  arma::field<arma::mat> out(10);
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
  	arma::mat V = gmwm_bootstrapper(GMWM.row(b), desc, max(tau), N, B, var_or_mu);
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
}
*/