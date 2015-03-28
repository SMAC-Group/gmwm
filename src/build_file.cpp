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

// [[Rcpp::export]]
unsigned int sum_field_vec(const arma::field<arma::vec>& x){
  unsigned int nelems = x.n_elem;
  unsigned int total_elems = 0;
  
  for(unsigned int i = 0; i < nelems; i++){
    total_elems += sum(x(i));
  }
  
  return total_elems;
}

// [[Rcpp::export]]
double mean_diff(const arma::vec& x){
  return arma::mean(diff_cpp(x, 1, 1));
}

// [[Rcpp::export]]
arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
    arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // AR 1
  	  if(desc[i] == "AR1" ){
        
        // Apply model specific parameter transformation
        if(model_type == "imu"){
  	      starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
        }else{ // O.W. SSM case
          starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
        }
        
        // Increase index for phi term
  	    ++i_theta;
        
        // Handle the SIGMA2 term
  	    starting(i_theta) = log(theta(i_theta));
        
  	  }
      else if( desc[i] == "ARMA" ) {
        
        arma::vec arma_params = objdesc(i);
        
        unsigned int p = arma_params(0); // AR(P)
        unsigned int q = arma_params(1); // MA(Q)
        
        unsigned int param_space = i_theta + p + q;
        
        // Change pseudo_logit to logit based on model type??
        starting.rows(i_theta, param_space - 1) = pseudo_logit(theta.rows(i_theta, param_space - 1));
        
        // Increment index
        i_theta += param_space;
        
        // Take care of SIGMA2 term
        starting(param_space) = log(theta(param_space));

      }
  	  else{
        starting(i_theta) = log(theta(i_theta));
  	  }
      ++i_theta;
  }  
    
  return starting;
}


// [[Rcpp::export]]
arma::colvec untransform_values(const arma::vec& theta, 
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
    arma::colvec result  = arma::zeros<arma::colvec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      
      std::string element_type = desc[i];
      // AR 1
  	  if(element_type == "AR1"){
        if(model_type == "imu"){
          result(i_theta) = arma::as_scalar(logit_inv(theta.row(i_theta)));
        }else{ // ssm
          result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
        }
  	    ++i_theta;
  	    result(i_theta) = exp(theta(i_theta));
  	  }
      else if(element_type == "ARMA" ) {
        
        arma::vec arma_params = objdesc(i);
        
        unsigned int p = arma_params(0); // AR(P)
        unsigned int q = arma_params(1); // MA(Q)
        
        unsigned int param_space = i_theta + p + q;
        
        // Change pseudo_logit_inv to logit_inv based on model type??
        result.rows(i_theta, param_space - 1) = pseudo_logit_inv(theta.rows(i_theta, param_space - 1));
        
        // Increment index to account for P+Q values
        i_theta += param_space;
        
        // Take care of SIGMA2 term
        result(param_space) = exp(theta(param_space));

      }
  	  else {
        result(i_theta) = exp(theta(i_theta));
  	  }
      
      ++i_theta;
  }  
    
  return result;
}


// hiding this function for the moment
double objFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec untransformed_theta = untransform_values(theta, desc, objdesc, model_type);
  arma::vec wv_theo = theoretical_wv(untransformed_theta, desc, objdesc, tau);
  arma::vec standardized = 1-wv_theo/wv_empir;
	// Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

double objFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega,const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec untransformed_theta = untransform_values(theta, desc, objdesc, model_type);

  arma::vec wv_theo = theoretical_wv(untransformed_theta, desc, objdesc, tau);

  // Compute quandratic form
	arma::vec dif = wv_theo - wv_empir;
	return arma::as_scalar(trans(dif)*omega*dif);
}

// [[Rcpp::export]]
double getObjFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

  return objFunStarting(transformed_theta, desc, objdesc, model_type, wv_empir, tau);
}

// [[Rcpp::export]]
double getObjFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega,const arma::vec& wv_empir, const arma::vec& tau){
  
    arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

    return objFun(transformed_theta, desc, objdesc, model_type, arma::inv(omega), wv_empir, tau);
}

/// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                          const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        _["desc"] = desc,
                        _["objdesc"] = objdesc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, 
                     const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                     const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun),
                        _["desc"] = desc,
                        _["objdesc"] = objdesc,
                        _["model_type"] = model_type,
                        _["omega"] = omega,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}



// [[Rcpp::export]]
arma::vec gmwm_bootstrapper(const arma::vec&  theta,
                            const std::vector<std::string>& desc, arma::field<arma::vec>& objdesc,
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

// [[Rcpp::export]]
arma::vec ar1_draw(unsigned int draw_id, double last_phi, double sigma_tot, std::string model_type){
  arma::vec temp(2);
  
  
  if(draw_id == 0){
    if(model_type == "imu"){
      // Draw from triangle distributions for phi
      double U = R::runif(0.0, 1.0/3.0);
      
      // Draw for phi
      temp(0) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
      temp(1) = R::runif(0.95*sigma_tot*(1-square(temp(0))), sigma_tot);
    }
    else{ // ssm
      // Draw for phi
      temp(0) = R::runif(-0.9999999999999, 0.9999999999999);
      // Draw for sigma
      temp(1) = R::runif(0.0000000000001, sigma_tot);
    }
  }
  else{
    
    if(draw_id!=1){
      // Draw for phi on i >= 3
      temp(0) = R::runif(last_phi,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    else{
      // Draw for phi on i==1
      temp(0) = R::runif(0.995,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    
    // Draw for process variance
    temp(1) = R::runif(0.0, 0.01*sigma_tot*(1-square(temp(0))) ); // VERIFY THIS SHOULD BE PHI VALUE!!
    
  } // end if
  
  return temp;
}

// [[Rcpp::export]]
unsigned int count_AR1s(std::vector<std::string> s) {
  unsigned int count = 0;

  for (unsigned int i = 0; i < s.size(); i++)
    if (s[i] == "AR1") count++;

  return count;
}

// [[Rcpp::export]]
std::map<std::string, int> count_models(const std::vector<std::string>& desc){  	
  std::map<std::string, int> w;	
  
  // We want to see the only the following objects with these initial values
  w["AR1"]=0;
  w["ARMA"]=0;
  w["DR"]=0;		
  w["RW"]=0;		
  w["QN"]=0;		
  w["WN"]=0;		
  
  for (unsigned int i = 0; i < desc.size(); i++) {		
    ++w[desc[i]];		
  }		
  
  return w;		
} 

//' @title Randomly guess a starting parameter
//' @description Sets starting parameters for each of the given parameters. 
//' @param desc A \code{vector<string>} that contains the model's components.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @param model_type A \code{string} that indicates whether it is an SSM or IMU.
//' @param num_params An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
//' @param expect_diff A \code{double} that contains the mean of the first difference of the data
//' @param N A \code{integer} that contains the number of observations in the data.
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
//' @param B A \code{integer} that indicates how many random draws that should be performed.
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec guess_initial(const std::vector<std::string>& desc, arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::vec& wv_empir, const arma::vec& tau, unsigned int B=1000){
                          
  // Obtain the sum of variances for sigma^2_total.
  double sigma_tot = arma::sum(wv_empir);
    
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
    
  double min_obj_value = std::numeric_limits<double>::max();
  
  std::map<std::string, int> models = count_models(desc);
  
  unsigned int num_desc = desc.size();
  
  unsigned int AR1_counter; // identifiability hack. =(
  double prev_phi; // ar1_draw needs external memory  
  
  // Generate parameters for the model
  for(unsigned int b = 0; b < B; b++){
    
    unsigned int i_theta = 0;

    if(models["WN"] >= 1 && model_type=="imu"){
      AR1_counter = 2;
      prev_phi = .9;
    }
    else{ // Multiple WN in model
      AR1_counter = 0;
      prev_phi = 0;
    }
    
        
    for(unsigned int i = 0; i < num_desc; i++){
      std::string element_type = desc[i];
      
      if(element_type == "AR1"){
        temp_theta.rows(i_theta, i_theta + 1) = ar1_draw(AR1_counter, prev_phi, sigma_tot, model_type);
        prev_phi = temp_theta(i_theta);
        i_theta++; // needed to account for two parameters (e.g. phi + sigma2). Second shift at end.
        AR1_counter++;
      }
      else if(element_type == "ARMA"){
        //  This needs to be implemented.
      }
      else if(element_type == "DR"){   
        if(expect_diff > 0){
          expect_diff = R::runif(0,2*expect_diff);
        }else{
          expect_diff = R::runif(2*expect_diff,0);
        }
        temp_theta(i_theta) = expect_diff;
      }
      else if(element_type == "QN"){
        temp_theta(i_theta) = R::runif(.0000001, sigma_tot);
      }
      else if(element_type == "RW"){
        temp_theta(i_theta) = R::runif(sigma_tot/double(N*1000.0), 2.0*sigma_tot/double(N));
      }
      else{ // WN
        temp_theta(i_theta) = R::runif(sigma_tot/2.0, sigma_tot);
      }
      i_theta ++;
    } // end for
    
    arma::vec tvalues = transform_values(temp_theta, desc, objdesc, model_type);
        
    double obj = objFunStarting(tvalues, desc, objdesc, model_type, wv_empir, tau);
    
    if(min_obj_value > obj){
      min_obj_value = obj;
      starting_theta = temp_theta;
    } //end if
  } // end for
  
  return starting_theta;
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
arma::rowvec gmwm_cpp(const arma::vec& theta,
                          const std::vector<std::string>& desc, arma::field<arma::vec>& objdesc, std::string model_type, 
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
                          const std::vector<std::string>& desc, arma::field<arma::vec>& objdesc, std::string model_type, 
                          const arma::mat& V, const arma::vec& wv_empir,
                          const arma::vec& tau){
                                 
  // Number of parameters
  //unsigned int num_param = theta.n_elem;
    
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