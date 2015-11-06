/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
 * (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#include <RcppArmadillo.h>
#include "process_to_wv.h"

// Required for ARMA to WV process
#include "rtoarmadillo.h"

// Need square
#include "inline_functions.h"

using namespace Rcpp;

/* ----------------------------- Start Process to WV Functions ------------------------------- */

//' @title ARMA process to WV
//' @description This function computes the (haar) WV of an ARMA process
//' @param ar A \code{vec} containing the coefficients of the AR process
//' @param ma A \code{vec} containing the coefficients of the MA process
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @param sigma A \code{double} containing the residual variance
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @examples
//' arma_to_wv(c(.23,.43), c(.34,.41,.59), 2^(1:9), 3)
//' @seealso \code{\link{ARMAtoMA_cpp}},\code{\link{ARMAacf_cpp}}
// [[Rcpp::export]]
arma::vec arma_to_wv(arma::vec ar, arma::vec ma, arma::vec tau, double sigma) {
  
  arma::vec n = arma::sort(tau/2);
  unsigned int ntau = tau.n_elem;
  double sig2 = (arma::sum(arma::square(ARMAtoMA_cpp(ar,ma,1000)))+1)*sigma;
  
  arma::vec acfvec = ARMAacf_cpp(ar, ma, as_scalar(tau.tail(1)) - 1);
  
  arma::vec wvar(ntau);
  
  unsigned int scale = n(0);
  
  // initial starting term
  double term4 = acfvec(scale);
  wvar(0)=( ( ( scale*(1.0-term4)) / square(scale))*sig2)/2.0;
  
  for(unsigned int j = 1; j < ntau; j++){
    scale = n(j);
    arma::vec boh(scale - 1);
    for (unsigned int i=1; i<= scale - 1; i++){
      double term1 = acfvec(scale-i);
      double term2 = acfvec(i);
      double term3 = acfvec(2*scale-i);
      // Account for starting loop at 1 instead of 0.
      boh(i-1)=i*((2.0*term1)-term2-term3);
    }
    double term4 = acfvec(scale);
    wvar(j)=((( (scale*(1.0-term4) ) + arma::sum(boh) ) /square(scale))*sig2)/2.0;
  }
  
  return wvar;
}

//' @title Helper Function for ARMA to WV Approximation
//' @description Indicates where the minimum ARMAacf value is and returns that as an index.
//' @param ar A \code{vec} containing the coefficients of the AR process
//' @param ma A \code{vec} containing the coefficients of the MA process
//' @param last_tau An \code{int} the Jth scale of 2^(1:J)
//' @param alpha A \code{double} indicating the cutoff.
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @keywords internal
//' @seealso \code{\link{arma_to_wv_app}}
// [[Rcpp::export]]
double acf_sum(arma::vec ar, arma::vec ma, unsigned int last_tau, double alpha = 0.99){
  arma::vec obj = abs(ARMAacf_cpp(ar,ma,last_tau -1));
  
  obj = abs(cumsum(obj)/sum(obj) - alpha);
  return as_scalar(find(obj == min(obj))) + 1;
}

//' @title ARMA process to WV approximation
//' @description This function computes the (haar) WV of an ARMA process
//' @param ar A \code{vec} containing the coefficients of the AR process
//' @param ma A \code{vec} containing the coefficients of the MA process
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @param sigma A \code{double} containing the residual variance
//' @param alpha A \code{double} indicating the cutoff.
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @keywords internal
//' @examples
//' arma_to_wv_app(c(.23,.43), c(.34,.41,.59), 2^(1:9), 3, .9)
//' @seealso \code{\link{ARMAtoMA_cpp}},\code{\link{ARMAacf_cpp}}
// [[Rcpp::export]]
arma::vec arma_to_wv_app(arma::vec ar, arma::vec ma, arma::vec tau, double sigma, double alpha = 0.9999) {
  
  arma::vec n = arma::sort(tau/2);
  unsigned int ntau = tau.n_elem;
  
  unsigned int lag = acf_sum(ar, ma, as_scalar(tau.tail(1)), alpha);
    
  double sig2 = (arma::sum(arma::square(ARMAtoMA_cpp(ar,ma,1000)))+1)*sigma;
  
  arma::vec acfvec = ARMAacf_cpp(ar, ma, lag);
  
  arma::vec wvar(ntau);
  
  
  
  double term1;
  double term2;
  double term3;
  
  unsigned int scale = n(0);
  
  // initial starting term
  double term4 = acfvec(scale);
  wvar(0)=( ( ( scale*(1.0-term4)) / square(scale))*sig2)/2.0;

  
  for(unsigned int j = 1; j < ntau; j++){
    scale = n(j);
    arma::vec boh(scale - 1);
    for (unsigned int i=1; i<= scale - 1; i++){
      if(scale-i > lag){
        term1 = 0;
      }else{
        term1 = acfvec(scale-i);
      }
      
      if(i > lag){
        term2 = 0;
      } else{
        term2 = acfvec(i);
      }
      
      if(2*scale - i > lag){
        term3 = 0;
      } else{
        term3 = acfvec(2*scale-i);
      }

      // Account for starting loop at 1 instead of 0.
      boh(i-1)=i*((2.0*term1)-term2-term3);
    }
    if(scale > lag){
      term4 = 0;
    }else{
      term4 = acfvec(scale);
    }
    wvar(j)=((( (scale*(1.0-term4) ) + arma::sum(boh) ) /square(scale))*sig2)/2.0;
  }
  
  return wvar;
}




//' @title Quantisation Noise to WV
//' @description This function compute the WV (haar) of a Quantisation Noise (QN) process
//' @param q2 A \code{double} corresponding to variance of drift
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the QN.
//' @examples
//' x.sim = 1:1000
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = qn_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec qn_to_wv(double q2, const arma::vec& tau){
  return 6.0*q2/arma::square(tau);
}

//' @title White Noise to WV
//' @description This function compute the WV (haar) of a White Noise process
//' @param sig2 A \code{double} corresponding to variance of WN
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the white noise.
//' @examples
//' x.sim = cumsum(rnorm(100000))
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = wn_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec wn_to_wv(double sig2, arma::vec tau){
  return sig2/tau;
}


//' @title Random Walk to WV
//' @description This function compute the WV (haar) of a Random Walk process
//' @param sig2 A \code{double} corresponding to variance of RW
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the random walk.
//' @examples
//' x.sim = cumsum(rnorm(100000))
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = rw_to_wv(1,tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec rw_to_wv(double sig2, const arma::vec& tau){
  return sig2*((2.0*arma::square(tau) + 4.0)/(24.0*tau));
}


//' @title Drift to WV
//' @description This function compute the WV (haar) of a Drift process
//' @param omega A \code{double} corresponding to variance of drift
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the drift.
//' @examples
//' x.sim = 1:1000
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = dr_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec dr_to_wv(double omega,const arma::vec& tau){
	return square(omega)*arma::square(tau)/16.0;
}

//' @title AR1 process to WV
//' @description This function compute the WV (haar) of an AR(1) process
//' @param phi A \code{double} that is the phi term of the AR(1) process
//' @param sig2 A \code{double} corresponding to variance of AR(1) process
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the AR(1) process.
//' @examples
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sigma2 = 4 )
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec ar1_to_wv(double phi, double sig2, const arma::vec& tau){
  unsigned int size_tau = tau.n_elem;
  arma::vec temp_term(size_tau);
  arma::vec temp_term_redux(size_tau);
  for(unsigned int i=0; i< size_tau; i++){
    temp_term(i) = 4*pow(phi,(tau(i)/2 + 1));
    temp_term_redux(i) = pow(phi,(tau(i)+1));
  }
	return ((tau/2.0 - 3.0*phi - tau/2.0*pow(phi,2) + temp_term - temp_term_redux)/(arma::square(tau/2.0)*pow(1-phi,2)*(1-pow(phi,2)))*sig2)/2.0;
}

//' @title Model Process to WV
//' @description This function computes the summation of all Processes to WV (haar) in a given model
//' @param theta A \code{vec} containing the list of estimated parameters.
//' @param desc A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @param tau A \code{vec} containing the scales e.g. 2^(1:J)
//' @return A \code{vec} containing the wavelet variance of the model.
//' @examples
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sigma2 = 4 )
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, tau)
//' plot(tau, wv.theo, col = "red")
//' @keywords internal
// [[Rcpp::export]]
arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau){
  
  unsigned int num_desc = desc.size();
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
    
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){
    
    // Add ARMA
  
    double theta_value = theta(i_theta);
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1"){

      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo += ar1_to_wv(theta_value, sig2, tau);
    }
    else if(element_type == "ARMA"){
      
      arma::vec model_params = objdesc(i);
      
      unsigned int p = model_params(0);
      unsigned int q = model_params(1);
      
      arma::vec ar;
      arma::vec ma;
      

      if(p == 0){
        ar = arma::zeros<arma::vec>(0);
      }else{
        ar = theta.rows(i_theta,i_theta+p-1);
      }
      
      i_theta += p;
      
      if(q == 0){
        ma = arma::zeros<arma::vec>(0); 
      }else{
        ma = theta.rows(i_theta,i_theta+q-1);
      }
      
      i_theta += q;

      
      wv_theo += arma_to_wv(ar, ma, tau, theta(i_theta));
    }
    // DR
    else if(element_type == "DR"){
      wv_theo += dr_to_wv(theta_value, tau);
    }
    // QN
    else if(element_type == "QN"){
      wv_theo += qn_to_wv(theta_value, tau);
    }
    // RW
    else if(element_type == "RW"){
      wv_theo += rw_to_wv(theta_value, tau);
    }
    // WN
    else{
      wv_theo += wn_to_wv(theta_value, tau);
    }
    
    ++i_theta;
  }

  return wv_theo;
}


//' @title Each Models Process Decomposed to WV
//' @description This function computes each process to WV (haar) in a given model.
//' @param theta A \code{vec} containing the list of estimated parameters.
//' @param desc A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @param tau A \code{vec} containing the scales e.g. 2^(1:J)
//' @return A \code{mat} containing the wavelet variance of each process in the model
//' @examples
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sigma2 = 4 )
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, tau)
//' plot(tau, wv.theo, col = "red")
//' @keywords internal
// [[Rcpp::export]]
arma::mat decomp_theoretical_wv(const arma::vec& theta, 
                                const std::vector<std::string>& desc,
                                const arma::field<arma::vec>& objdesc, const arma::vec& tau){
  
  unsigned int num_desc = desc.size();
  arma::mat wv_theo = arma::zeros<arma::mat>(tau.n_elem, num_desc);
    
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){
    
    // Add ARMA
  
    double theta_value = theta(i_theta);
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1"){

      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo.col(i) = ar1_to_wv(theta_value, sig2, tau);
    }
    else if(element_type == "ARMA"){
      
      arma::vec model_params = objdesc(i);
      
      unsigned int p = model_params(0);
      unsigned int q = model_params(1);
      
      arma::vec ar;
      arma::vec ma;
      

      if(p == 0){
        ar = arma::zeros<arma::vec>(0);
      }else{
        ar = theta.rows(i_theta,i_theta+p-1);
      }
      
      i_theta += p;
      
      if(q == 0){
        ma = arma::zeros<arma::vec>(0); 
      }else{
        ma = theta.rows(i_theta,i_theta+q-1);
      }
      
      i_theta += q;

      
      wv_theo.col(i) = arma_to_wv(ar, ma, tau, theta(i_theta));
    }
    // DR
    else if(element_type == "DR"){
      wv_theo.col(i) = dr_to_wv(theta_value, tau);
    }
    // QN
    else if(element_type == "QN"){
      wv_theo.col(i) = qn_to_wv(theta_value, tau);
    }
    // RW
    else if(element_type == "RW"){
      wv_theo.col(i) = rw_to_wv(theta_value, tau);
    }
    // WN
    else{
      wv_theo.col(i) = wn_to_wv(theta_value, tau);
    }
    
    ++i_theta;
  }

  return wv_theo;
}

//' @title Decomposed WV to Single WV
//' @description This function computes the combined processes to WV (haar) in a given model.
//' @param decomp A \code{mat} with scales as rows and processes as columns
//' @return A \code{vec} containing the wavelet variance of the process for the overall model
//' @examples
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sigma2 = 4 )
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, tau)
//' plot(tau, wv.theo, col = "red")
//' @keywords internal
// [[Rcpp::export]]
arma::vec decomp_to_theo_wv(const arma::mat& decomp){

  return arma::sum(decomp, 1);
}

/* ------------------------------ End Process to WV Functions --------------------------------- */
