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

// Inference goodies
#include "inference.h"

// Model Selection criterion
#include "analytical_matrix_derivatives.h"

// Model Selection criterion
#include "model_selection.h"

// Count Models
#include "ts_checks.h"

// arima
#include "arima_gmwm.h"

using namespace arma;
using namespace Rcpp;

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
                            unsigned int H, bool diagonal_matrix){
  unsigned int nb_level = floor(log2(N));
    
  arma::mat res(H, nb_level);
	for(unsigned int i=0; i<H; i++){
	  
	  // Generate x_t ~ F_theta
		arma::vec x = gen_model(N, theta, desc, objdesc);

    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level, "periodic");
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter(), "modwt");
  
    // Obtain WV
		arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
    // Store WV into matrix
	  res.row(i) = arma::trans(wv_x);
	}
  
  // Do we need a diagnoal covariance matrix? 
  if(diagonal_matrix){
    return arma::diagmat(arma::cov(res));
  }
  
  // Return the full covariance matrix
	return arma::cov(res);
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
//'   \item{"ARMA"}{an autoregressiveß moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
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
//' @examples
//' # Coming soon
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
arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, unsigned int N, double expect_diff, 
                                      arma::mat orgV, arma::vec scales, arma::vec wv_empir,
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
  
  // Do we need to run a guessing algorithm?
  if(starting){

    theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wv_empir, scales, G);
    
    guessed_theta = theta;
  }

  // Obtain the GMWM estimator estimates.
  theta = gmwm_engine(theta, desc, objdesc, model_type, 
                      wv_empir, omega, scales, starting);
  
  // Bootstrap the V matrix
  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        V = gmwm_bootstrapper(theta, desc, objdesc, N, robust, eff, H);
        omega = arma::inv(diagmat(V));
        theta = gmwm_engine(theta, desc, objdesc, model_type,
                      wv_empir, omega, scales, starting);
    }
  }

  // Obtain the theoretical WV.
  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);

  // Export calculations to R.
  arma::field<arma::mat> out(5);
  out(0) = theta;
  out(1) = guessed_theta;
  out(2) = V;
  out(3) = theo;
  out(4) = decomp_theo;
  
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
//' @param inference A \code{bool} that indicates whether inference should be run on the supplied model.
//' @return A \code{field<mat>} that contains a list of ever-changing estimates...
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
arma::field<arma::mat> gmwm_master_cpp(const arma::vec& data, 
                                      arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, bool starting,
                                      double alpha, 
                                      std::string compute_v, unsigned int K, unsigned int H,
                                      unsigned int G, 
                                      bool robust, double eff, bool inference){
  // Variable Declarations
  
  // Length of the Time Series
  unsigned int N = data.n_elem;
  
  // Number of Scales (J)
  unsigned int nlevels = floor(log2(N));
  
  // Number of parameters
  unsigned int np = theta.n_elem;
  
  // Obtain the diagonal matrix when bootstrapping
  bool diagonal_matrix = true;
  
  // Take the mean of the first difference
  double expect_diff = mean_diff(data);
  
  // Guessed values of Theta (user supplied or generated)
  arma::vec guessed_theta = theta;
  
  // MODWT decomp
  arma::field<arma::vec> modwt_decomp = modwt_cpp(data, "haar", nlevels, "periodic");
  
  // Obtain WV and confidence intervals
  arma::mat wvar = wvar_cpp(modwt_decomp, robust, eff, alpha, "eta3", "haar");
  
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
  
  // Guess starting values for the theta parameters
  if(starting){
    if(desc[0] == "ARMA" && desc.size() == 1){
      
      theta = Rcpp_ARIMA(data, objdesc(0)); 
      starting = false;
    }else{     
      theta = guess_initial(desc, objdesc, model_type, np, expect_diff, N, wv_empir, scales, G);
    }
    guessed_theta = theta;
  }

  Rcpp::Rcout << "Guessed" << guessed_theta << std::endl;
  // Obtain the GMWM estimator's estimates.
  theta = gmwm_engine(theta, desc, objdesc, model_type, 
                      wv_empir, omega, scales, starting);


  // Set values to compute inference
  if(inference){
    diagonal_matrix = false;
    compute_v = "bootstrap";
  }

  // Enable bootstrapping
  if(compute_v == "bootstrap"){
    for(unsigned int k = 0; k < K; k++){
        V = gmwm_bootstrapper(theta, desc, objdesc, N, robust, eff, H, diagonal_matrix);
        omega = arma::inv(diagmat(V));
        theta = gmwm_engine(theta, desc, objdesc, model_type, wv_empir, omega, scales, starting);
    }
  }
  
  // Obtain the objective value function
  arma::vec obj_value(1);
  obj_value(0) = getObjFun(theta, desc, objdesc,  model_type, omega, wv_empir, scales); 
  
  // Decomposition of the WV.
  arma::mat decomp_theo = decomp_theoretical_wv(theta, desc, objdesc, scales);
  arma::vec theo = decomp_to_theo_wv(decomp_theo);

  
  // Inference test result storage
  arma::vec gof_test;
  arma::mat ci_inf; 
  arma::vec score;

  // Generate inference information
  if(inference){
    
    // Take derivatives
    arma::mat D = derivative_first_matrix(theta, desc, objdesc, scales);
    arma::mat At_j = derivative_second_matrix(theta, desc, objdesc, scales);
    
    // Obtain a confidence interval for the parameter estimates AND calculate chisq goodness of fit
    arma::field<arma::mat> cat = inference_summary(theta, desc,  objdesc, model_type, scales,
                                                   D, V, omega, wv_empir, alpha);

    // Obtain the difference
    arma::vec diff = theo - wv_empir;
    
    // Calculate the model score according to model selection criteria paper
    score = model_score(D, At_j, omega, V,  diff, N);
    
    
    ci_inf = cat(0);
    gof_test = cat(1);
    
  }
  
  // Export information back
  arma::field<arma::mat> out(15);
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
  out(10) = ci_inf;
  out(11) = gof_test;
  out(12) = obj_value;
  out(13) = score;
  out(14) = omega;
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