#include <RcppArmadillo.h>
#include <string>
#include <map>
// #include <omp.h>

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

using namespace arma;
using namespace Rcpp;

// // [[Rcpp::plugins(openmp)]]

// computes theoretical wv
inline arma::vec theoretical_wv(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                                const arma::vec& wv_empir,
                                const arma::vec& tau, int N){
  
  unsigned int num_desc = desc.size();
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
    
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){
    
    // Add ARMA
  
    double theta_transformed = exp(theta(i_theta));

    // AR 1
    if(desc[i] == "AR1"){
      if(model_type == "imu"){
        theta_transformed = arma::as_scalar(logit_inv(theta.row(i_theta)));
      }
      else{ // ssm
        theta_transformed = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
      }
      ++i_theta;
      double sig2 = exp(theta(i_theta));
      
      // Compute theoretical WV
      wv_theo += ar1_to_wv(theta_transformed, sig2, tau);
    }
    // RW
    else if(desc[i] == "RW"){
      wv_theo += rw_to_wv(theta_transformed, tau);
    }
    // DR
    else if(desc[i] == "DR"){
      wv_theo += dr_to_wv(theta_transformed, tau);
    }
    // WN
    else{
      wv_theo += wn_to_wv(theta_transformed, tau);
    }
    ++i_theta;
  }

  return wv_theo;
}

// hiding this function for the moment
double objFunStarting(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, model_type, wv_empir, tau, N);
  arma::vec standardized = 1-wv_theo/wv_empir;
	// Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

double objFun(const arma::vec& theta, const arma::mat& omega,
                 const std::vector<std::string>& desc, std::string model_type,
                 const arma::vec& wv_empir, const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, model_type, wv_empir, tau, N);

  // Compute quandratic form
	arma::vec dif = wv_theo - wv_empir;
	return arma::as_scalar(trans(dif)*omega*dif);
}

/// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta, const std::vector<std::string>& desc, std::string model_type,
                          const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        _["desc"] = desc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, const std::vector<std::string>& desc, std::string model_type,
                     const arma::mat& omega, const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun),
                        _["omega"] = omega,
                        _["desc"] = desc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}


//////////// ARMA SECTION


// hiding this function for the moment
double objFunStarting_ARMA(const arma::vec& theta, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
  arma::vec ar;
  arma::vec ma;

  if(p == 0){
    ar = arma::zeros<arma::vec>(0);
  }else{
    ar = theta.rows(0,p-1);
  }
  
  if(q == 0){
    ma = arma::zeros<arma::vec>(0); 
  }else{
    ma = theta.rows(p,p+q-1);
  }
  
  double sigma2 = theta(theta.n_elem-1);
  arma::vec wv_theo = arma_to_wv(ar, ma, tau, sigma2);

  arma::vec standardized = 1-wv_theo/wv_empir;
  // Compute quandratic form
  return arma::as_scalar(trans(standardized)*(standardized));
}

arma::vec Rcpp_OptimStart_ARMA(const arma::vec& theta, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting_ARMA),
                        _["p"] = p,
                        _["q"] = q,
                        _["tau"] = tau,
                        _["wv_empir"] = wv_empir);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

double objFun_ARMA(const arma::vec& theta, const arma::mat& omega,  int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
                   
  arma::vec ar;
  arma::vec ma;
  
  if(p == 0){
    ar = arma::zeros<arma::vec>(0);
  }else{
    ar = pseudo_logit_inv(theta.rows(0,p-1));
  }
  
  if(q == 0){
    ma = arma::zeros<arma::vec>(0); 
  }else{
    ma = pseudo_logit_inv(theta.rows(p,p+q-1));
  }
  
  arma::vec wv_theo = arma_to_wv(ar, ma, tau, exp(theta(theta.n_elem-1)) );
  // Compute quandratic form
  arma::vec dif = wv_theo - wv_empir;
  return arma::as_scalar(trans(dif)*omega*dif);
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim_ARMA(const arma::vec& theta, const arma::mat& omega, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun_ARMA),
                        _["omega"] = omega,
                        _["p"] = p,
                        _["q"] = q,
                        _["tau"] = tau,
                        _["wv_empir"] = wv_empir);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}


///// END ARMA


/// [[Rcpp::export]]
arma::vec set_starting_values(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type){
    arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
  	  // AR 1
  	  if(desc[i] == "AR1" ){
        if(model_type == "imu"){
  	      starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
        }else{ // ssm
          starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
        }
  	    ++i_theta;
  	    starting(i_theta) = log(theta(i_theta));
  	    ++i_theta;
  	  }
  	  else{
        starting(i_theta) = log(theta(i_theta));
  	    ++i_theta;
  	  }
  }  
    
  return starting;
}


// [[Rcpp::export]]
arma::rowvec set_result_values(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type){
    arma::rowvec result  = arma::zeros<arma::rowvec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // AR 1
  	  if(desc[i] == "AR1"){
        if(model_type == "imu"){
          result(i_theta) = arma::as_scalar(logit_inv(theta.row(i_theta)));
        }else{ // ssm
          result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
        }
  	    ++i_theta;
  	    result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
  	  else {
        result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
  }  
    
  return result;
}


/// [[Rcpp::export]]
arma::vec set_starting_values_arma(arma::vec theta, int p, int q){
  theta.rows(0,p+q-1) = pseudo_logit(theta.rows(0,p+q-1));
  theta(theta.n_elem -1) = log(theta(theta.n_elem -1));
  
  return theta;
}


/// [[Rcpp::export]]
arma::vec set_result_values_arma(arma::vec theta, int p, int q){    
  theta.rows(0,p+q-1) = pseudo_logit_inv(theta.rows(0,p+q-1));
  theta(theta.n_elem -1) = exp(theta(theta.n_elem -1));
  
  return theta;
}


// [[Rcpp::export]]
arma::vec gmwm_bootstrapper(const arma::vec&  theta, const std::vector<std::string>& desc, 
                            unsigned int tau, unsigned int N, bool robust, double eff,
                            unsigned int B = 100){
  unsigned int nb_level = floor(log2(N));
  	
	arma::mat res(B, tau+1);
	for(unsigned int i=0; i<B; i++){
		arma::vec x = gen_model(N, theta, desc);

    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level, "periodic");
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter(), "modwt");
  
		arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
	  res.row(i) = arma::trans(wv_x);
	}
	return cov(res);
}

inline arma::vec ar1_draw(unsigned int num_ars, double sigma_tot, std::string model_type){
  
  unsigned int num_params = 2*num_ars;
  arma::vec temp(num_params);
  
  for(unsigned int i = 0; i < num_ars; i++){
    
    // Draw from triangle distributions for phi
    double U = R::runif(0.0, 1.0/3.0);
    
    if(i == 0){
      if(model_type == "imu"){
        // Draw for phi
        temp(2*i) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
        temp(2*i+1) = R::runif(0.95*sigma_tot*(1-square(temp(2*i))), sigma_tot);
      }
      else{ // ssm
        // Draw for phi
        temp(2*i) = R::runif(-0.9999999999999, 0.9999999999999);
        // Draw for sigma
        temp(2*i+1) = R::runif(0.0000000000001, sigma_tot);
      }
    }
    else{
      
      if(i!=1){
          // Draw for phi on i >= 3
          temp(2*i) = R::runif(temp(2*(i-1)),0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
      }
      else{
          // Draw for phi on i==1
          temp(2*i) = R::runif(0.995,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
      }
      
      // Draw for process variance
      temp(2*i+1) = R::runif(0.0, 0.01*sigma_tot*(1-square(temp(2*i+1))) );
    } // end if
    
  } // end for
  
  return temp;
}

inline arma::vec unif_sigma_sample(unsigned int num, double start, double end){
  arma::vec temp(num);
  
  for(unsigned int i = 0; i<num; i++){
    temp(i) = R::runif(start,end);
  }
  
  return temp;
}


// @title Randomly guess a starting parameter
// @description Sets starting parameters for each of the given parameters. 
// @usage guess_initial(signal, w, desc, model_type, num_param, wv_empir, tau, N, B)
// @param signal A \code{vec} that contains the data
// @param w A \code{map<string,int>} that lists supported models and the amount in the model.
// @param model_type A \code{string} that indicates whether it is an SSM or IMU.
// @param num_params An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
// @param wv_empir A \code{vec} that contains the empirical wavelet variance.
// @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
// @param N A \code{integer} that indicates the signal length
// @param B A \code{integer} that indicates how many random draws that should be performed.
// @return A \code{vec} containing smart parameter starting guesses to be iterated over.
// @name guess_initial
// @docType methods
// @rdname guess_initial-methods
arma::vec guess_initial(arma::vec signal, std::map< std::string ,int>& w, const std::vector<std::string>& desc, std::string model_type,
                        unsigned int num_param, const arma::vec& wv_empir, const arma::vec& tau, 
                        unsigned int N, unsigned int B=1000){
  
  // Obtain the sum of variances for sigma^2_total.
  double sigma_tot = arma::sum(wv_empir);
    
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
    
  double min_obj_value = std::numeric_limits<double>::max();
  
  // Generate parameters for the model
  for(unsigned int b = 0; b < B; b++){
    
    unsigned int i_theta = 0;
    
    for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
      int out = p->second;
      if(out > 0){
        std::string type = p->first;
        if(type == "AR1"){
          temp_theta.rows(i_theta, i_theta + 2*out - 1) = ar1_draw(out, sigma_tot, model_type);
          i_theta += 2*out;
        }
        else if(type == "DR"){   
          double dr_ed = mean(diff_cpp(signal));
          if(dr_ed > 0){
            dr_ed = R::runif(0,2*dr_ed);
          }else{
            dr_ed = R::runif(2*dr_ed,0);
          }
          temp_theta.rows(i_theta, i_theta + out - 1).fill( dr_ed  );
          i_theta += out;
        }
        else if(type == "QN"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, .0000001, sigma_tot);
          i_theta += out;
        }
        else if(type == "RW"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot/(signal.n_elem*1000.0), 2.0*sigma_tot/signal.n_elem);
          i_theta += out;
        }
        else{ // WN
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot/2.0, sigma_tot);
          i_theta += out;
        }
      } // end if
    } // end for

    double obj = objFunStarting(temp_theta, desc, model_type, wv_empir, tau, N);
    
    if(min_obj_value > obj){
      min_obj_value = obj;
      starting_theta = temp_theta;
    } //end if
  } // end for
  
  return starting_theta;
}

// Counts 
inline std::map<std::string, int> counted_map(const std::vector<std::string>& desc){
  std::map<std::string, int> w;
  w["AR1"]=0;
  w["DR"]=0;
  w["RW"]=0;
  w["QN"]=0;
  w["WN"]=0;

  for (unsigned int i = 0; i < desc.size(); i++) {
        ++w[desc[i]];
  }
  
  return w;
} 

// 
inline unsigned int count_params(std::map<std::string, int>& w) {
  unsigned int num_params = 0;     
  for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
    std::string type = p->first;
    int num_models = p->second;
    if(type != "AR1" && num_models != 0){
      ++num_params;
    }
    else{
      if(num_models > 0){
        num_params += (2*num_models);
      }
    }
  }
  return num_params;
}

//' @title Count Number of Models (Alphanumeric)
//' @description Return a model count
//' @usage count_models_alpha(desc)
//' @param desc A \code{vector<string>} that contains the model type
//' @return A \code{vec} with the model counts.
//' @details 
//' The types of models supported are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' @examples
//' count_models_alpha(c("AR1","DR","AR1","WN"))
// [[Rcpp::export]]
arma::vec count_models_alpha(const std::vector<std::string>& desc) {
  
    std::map<std::string, int> w = counted_map(desc);

    arma::vec num_models = arma::zeros<arma::vec>(w.size());
    
    for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
              int out = p->second;
              num_models(std::distance(w.begin(), p)) = out;
    }
    return num_models;
}

// [[Rcpp::export]]
unsigned int num_model_params(const std::vector<std::string>& desc) {
  // Count the number of models we are working with
  std::map<std::string, int> w = counted_map(desc);
  
  // Return the total number of parameters we need to setup.
  return count_params(w);
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
arma::rowvec adv_gmwm_imu_ssm_cpp(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                   const arma::mat& V, const arma::vec& wv_empir,
                   const arma::vec& tau, unsigned int N){
                                 
  // Number of parameters
  //unsigned int num_param = theta.n_elem;
      
  // Starting values
  arma::vec starting_theta = set_starting_values(theta, desc, model_type);
  
  // Optimize Starting values via Jannick's Method
  starting_theta = Rcpp_OptimStart(starting_theta, desc, model_type, tau, wv_empir, N);

  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
      
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, model_type, omega, tau, wv_empir, N);
  
  return set_result_values(estim_GMWM, desc, model_type);                          
}


//' @title GMWM for IMU and SSM
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model.
//' @usage gmwm_imu_ssm_cpp(desc, signal, model_type, V, wv_empir, tau, N, B = 1000)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param model_type A \code{character string} indicating if the function should estimate an ARMA model ("ARMA"), a model for IMU sensor calibration ("IMU") or a state-space model ("SSM")
//' @param params A \code{vector} being numeric (if type = "ARMA") or character string (if type = "IMU" or type = "SSM")
//' @param robust A \code{bool} indicating if the function should provide a robust estimation of the model parameters (by default = FALSE).
//' @return gmwm A \code{list} that contains:
//' \itemize{
//'  \item{par}{The estimated model parameters}
//'  \item{CI}{The 95\% confidence intervals for the estimated model parameters.}
//' }
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). If type = "IMU" or "SSM", then
//' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
//' options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
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
arma::rowvec gmwm_imu_ssm_cpp(const std::vector<std::string>& desc, const arma::vec& signal, std::string model_type,
               const arma::mat& V, const arma::vec& wv_empir,
               const arma::vec& tau, unsigned int N, unsigned int B = 1000){
  
  // Count the number of models we are working with
  std::map<std::string, int> w = counted_map(desc);
  
  // Return the total number of parameters we need to setup.
  unsigned int num_param = count_params(w);
  
  // Give it a guess
  arma::vec guess_me = guess_initial(signal, w, desc, model_type, num_param, wv_empir, tau, N, B);

  // And return value...
  return adv_gmwm_imu_ssm_cpp(guess_me, desc, model_type, V, wv_empir, tau, N);                        
}

//' @title GMWM ARMA
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model based on ARMA.
//' @usage gmwm_arma_cpp(theta, V, p, q, tau, wv_empir)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param type A \code{character string} indicating if the function should estimate an ARMA model ("ARMA"), a model for IMU sensor calibration ("IMU") or a state-space model ("SSM")
//' @param params A \code{vector} being numeric (if type = "ARMA") or character string (if type = "IMU" or type = "SSM")
//' @param robust A \code{bool} indicating if the function should provide a robust estimation of the model parameters (by default = FALSE).
//' @return gmwm A \code{list} that contains:
//' \itemize{
//'  \item{par}{The estimated model parameters}
//'  \item{CI}{The 95\% confidence intervals for the estimated model parameters.}
//' }
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). 
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec gmwm_arma_cpp(const arma::vec& theta, const arma::mat& V, unsigned int p, unsigned int q,
                const arma::vec& tau, const arma::vec& wv_empir){
    
  arma::vec starting_theta = set_starting_values_arma(theta, p, q);

  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  arma::vec values = Rcpp_Optim_ARMA(theta, omega, p, q, tau, wv_empir);
  
  // Find GMWM estimator
    // Initialize it
  arma::rowvec GMWM = arma::trans(set_result_values_arma(values, p, q));
      
  return GMWM;                          
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