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
#include "process_to_wv.h"

// Required for ARMA to WV process
#include "rtoarmadillo.h"

// Need square
#include "inline_functions.h"

// Needed for sarma model support
#include "sarma.h"

/* ----------------------------- Start Process to WV Functions ------------------------------- */

//' ARMA process to WV
//' 
//' This function computes the Haar Wavelet Variance of an ARMA process
//' @param ar     A \code{vec} containing the coefficients of the AR process
//' @param ma     A \code{vec} containing the coefficients of the MA process
//' @param sigma2 A \code{double} containing the residual variance
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @details
//' The function is a generic implementation that requires a stationary theoretical autocorrelation function (ACF)
//' and the ability to transform an ARMA(\eqn{p},\eqn{q}) process into an MA(\eqn{\infty}{infinity}) (e.g. infinite MA process).
//' @template to_wv/haar_arma
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @examples
//' # Calculates the Haar WV for an ARMA(2,3).
//' wv.theo = arma_to_wv(c(.23,.43), c(.34,.41,.59), 3, 2^(1:9))
//' @seealso \code{\link{ARMAtoMA_cpp}}, \code{\link{ARMAacf_cpp}}, and \code{\link{arma11_to_wv}}
// [[Rcpp::export]]
arma::vec arma_to_wv(arma::vec ar, arma::vec ma, double sigma2, arma::vec tau) {
  
  arma::vec n = arma::sort(tau/2);
  unsigned int ntau = tau.n_elem;
  double sig2 = (arma::sum(arma::square(ARMAtoMA_cpp(ar,ma,1000)))+1)*sigma2;
  
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

//' ARMA process to WV Approximation
//' 
//' This function computes the (haar) WV of an ARMA process
//' @param ar A \code{vec} containing the coefficients of the AR process
//' @param ma A \code{vec} containing the coefficients of the MA process
//' @param sigma2 A \code{double} containing the residual variance
//' @template misc/tau
//' @param alpha A \code{double} indicating the cutoff.
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @keywords internal
//' @details
//' This function provides an approximation to the \code{\link{arma_to_wv}} as computation times
//' were previously a concern. However, this is no longer the case and, thus, this has been left
//' in for the curious soul to discover... 
//' @template to_wv/haar_arma
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @examples
//' # Performs an approximation of the Haar WV for an ARMA(2,3).
//' wv.theo = arma_to_wv_app(c(.23,.43), c(.34,.41,.59), 3, 2^(1:9), .9)
//' @seealso \code{\link{ARMAtoMA_cpp}}, \code{\link{ARMAacf_cpp}}, \code{\link{acf_sum}} and \code{\link{arma_to_wv}}
// [[Rcpp::export]]
arma::vec arma_to_wv_app(arma::vec ar, arma::vec ma, double sigma2, arma::vec tau, double alpha = 0.9999) {
  
  arma::vec n = arma::sort(tau/2);
  unsigned int ntau = tau.n_elem;
  
  unsigned int lag = acf_sum(ar, ma, as_scalar(tau.tail(1)), alpha);
    
  double sig2 = (arma::sum(arma::square(ARMAtoMA_cpp(ar,ma,1000)))+1)*sigma2;
  
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


//' ARMA(1,1) to WV
//' 
//' This function computes the WV (haar) of an Autoregressive Order 1 - Moving Average Order 1 (ARMA(1,1)) process.
//' @param phi    A \code{double} corresponding to the autoregressive term.
//' @param theta  A \code{double} corresponding to the moving average term. 
//' @param sigma2 A \code{double} the variance of the process. 
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the ARMA(1,1) process.
//' @details 
//' This function is significantly faster than its generalized counter part
//' \code{\link{arma_to_wv}}
//' 
//' @template to_wv/haar_arma11
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @seealso \code{\link{arma_to_wv}}
//' @examples
//' ntau = 7
//' tau = 2^(1:ntau)
//' wv.theo = arma11_to_wv(0.3, 0.1, 1, tau)
// [[Rcpp::export]]
arma::vec arma11_to_wv(double phi, double theta, double sigma2, const arma::vec& tau){
  
  unsigned int size_tau = tau.n_elem;
  arma::vec phi_tau_ov2(size_tau);
  arma::vec phi_tau(size_tau);
  for(unsigned int i=0; i< size_tau; i++){
    phi_tau_ov2(i) = pow(phi,tau(i)/2.0);
    phi_tau(i) = pow(phi,tau(i));
  }

  return (-2.0*sigma2*((-(theta + phi))*(1.0 + theta*phi)*
          (3.0 - 4.0*phi_tau_ov2 + phi_tau) - 0.5*square(1.0 + theta)*(-1.0 + square(phi))*tau)) /
          (std::pow(-1.0 + phi,3.0)*(1.0 + phi)*arma::square(tau));
}


//' AR(1) process to WV
//' 
//' This function computes the Haar WV of an AR(1) process
//' @param phi    A \code{double} that is the phi term of the AR(1) process
//' @param sigma2 A \code{double} corresponding to variance of AR(1) process
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the AR(1) process.
//' @details 
//' This function is significantly faster than its generalized counter part
//' \code{\link{arma_to_wv}}.
//' 
//' @template to_wv/haar_ar1
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @seealso \code{\link{arma_to_wv}}, \code{\link{arma11_to_wv}}
//' @examples
//' ntau = 7
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(.63, 1, tau)
// [[Rcpp::export]]
arma::vec ar1_to_wv(double phi, double sigma2, const arma::vec& tau){
  unsigned int size_tau = tau.n_elem;
  arma::vec temp_term(size_tau);
  arma::vec temp_term_redux(size_tau);
  for(unsigned int i=0; i< size_tau; i++){
    temp_term(i) = 4*pow(phi,(tau(i)/2 + 1));
    temp_term_redux(i) = pow(phi,(tau(i)+1));
  }
  return ((tau/2.0 - 3.0*phi - tau/2.0*pow(phi,2) + temp_term - temp_term_redux)/(arma::square(tau/2.0)*pow(1-phi,2)*(1-pow(phi,2)))*sigma2)/2.0;
}


//' Moving Average Order 1 (MA(1)) to WV
//' 
//' This function computes the WV (haar) of a Moving Average order 1 (MA1) process.
//' @param theta A \code{double} corresponding to the moving average term. 
//' @param sigma2  A \code{double} the variance of the process. 
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the MA(1) process.
//' @details 
//' This function is significantly faster than its generalized counter part
//' \code{\link{arma_to_wv}}.
//' 
//' @template to_wv/haar_ma1
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @seealso \code{\link{arma_to_wv}}, \code{\link{arma11_to_wv}}
//' @examples
//' ntau = 7
//' tau = 2^(1:ntau)
//' wv.theo = ma1_to_wv(.3, 1, tau)
// [[Rcpp::export]]
arma::vec ma1_to_wv(double theta, double sigma2, const arma::vec& tau){
  return sigma2 * (square(theta + 1.0) * tau - 6.0 * theta)/arma::square(tau);
}

//' Quantisation Noise (QN) to WV
//' 
//' This function compute the Haar WV of a Quantisation Noise (QN) process
//' @param q2  A \code{double} corresponding to variance of drift
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the QN.
//' @template to_wv/haar_qn
//' @template misc/haar_wv_formulae_link
//' @backref src/process_to_wv.cpp
//' @backref src/process_to_wv.h
//' @examples
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = qn_to_wv(.42, tau)
// [[Rcpp::export]]
arma::vec qn_to_wv(double q2, const arma::vec& tau){
  return 6.0*q2/arma::square(tau);
}

//' @title Gaussian White Noise to WV
//' @description This function compute the Haar WV of a Gaussian White Noise process
//' @param sigma2 A \code{double} corresponding to variance of WN
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the white noise.
//' @template to_wv/haar_wn
//' @template misc/haar_wv_formulae_link
//' @examples
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = wn_to_wv(1, tau)
// [[Rcpp::export]]
arma::vec wn_to_wv(double sigma2, arma::vec tau){
  return sigma2/tau;
}


//' @title Random Walk to WV
//' @description This function compute the WV (haar) of a Random Walk process
//' @param gamma2 A \code{double} corresponding to variance of RW
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the random walk.
//' @template to_wv/haar_rw
//' @template misc/haar_wv_formulae_link
//' @examples
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = rw_to_wv(.37, tau)
// [[Rcpp::export]]
arma::vec rw_to_wv(double gamma2, const arma::vec& tau){
  return gamma2*((arma::square(tau) + 2.0)/(12.0*tau));
}


//' @title Drift to WV
//' @description This function compute the WV (haar) of a Drift process
//' @param omega A \code{double} corresponding to the slope of the drift
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the drift.
//' @template to_wv/haar_dr
//' @template misc/haar_wv_formulae_link
//' @examples
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = dr_to_wv(-2.3, tau)
// [[Rcpp::export]]
arma::vec dr_to_wv(double omega, const arma::vec& tau){
	return square(omega)*arma::square(tau)/16.0;
}

//' Model Process to WV
//' 
//' This function computes the summation of all Processes to WV (haar) in a given model
//' @param theta   A \code{vec} containing the list of estimated parameters.
//' @param desc    A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @template misc/tau
//' @return A \code{vec} containing the wavelet variance of the model.
//' @template misc/haar_wv_formulae_link
//' @examples
//' model = AR1(.3,2) + RW(.21) + DR(.001)
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = theoretical_wv(model$theta, model$desc, model$objdesc, tau)
//' @keywords internal
// [[Rcpp::export]]
arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau){
  
  unsigned int num_desc = desc.size();
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
    
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){
  
    double theta_value = theta(i_theta);
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1" || element_type == "GM"){

      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo += ar1_to_wv(theta_value, sig2, tau);
    }
    else if(element_type == "MA1"){
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo += ma1_to_wv(theta_value, sig2, tau);
    }
    // WN
    else if(element_type == "WN"){
      wv_theo += wn_to_wv(theta_value, tau);
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
    // ARMA11
    else if(element_type == "ARMA11"){
      ++i_theta;
      double th = theta(i_theta);
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      wv_theo += arma11_to_wv(theta_value, th, sig2, tau);
    }
    else { // ARMA
      
      // Unpackage ARMA model parameter
      arma::vec model_params = objdesc(i);
      
      unsigned int pop = arma::sum(model_params.rows(0,3));
      
      // Extract theta_values (includes the active theta_value)
      // Takes np + nq + nsp + nsq values
      arma::vec theta_values = theta.rows(i_theta, i_theta + pop - 1);
      
      // Increase the theta counter
      i_theta += pop;
      
      // Setup parameters
      arma::field<arma::vec> psetup = sarma_expand(theta_values, model_params);
      
      // Pip into the gen_arima function!
      // Note this floors the function at d. 
      
      wv_theo += arma_to_wv(psetup(0), psetup(1), theta(i_theta), tau);
    }
    
    ++i_theta;
  }

  return wv_theo;
}


//' Each Models Process Decomposed to WV
//' 
//' This function computes each process to WV (haar) in a given model.
//' @param theta   A \code{vec} containing the list of estimated parameters.
//' @param desc    A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @template misc/tau
//' @return A \code{mat} containing the wavelet variance of each process in the model
//' @template misc/haar_wv_formulae_link
//' @examples
//' model = AR1(.3,2) + DR(.001)
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = decomp_theoretical_wv(model$theta, model$desc, model$objdesc, tau)
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
    if(element_type == "AR1" || element_type == "GM"){

      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo.col(i) = ar1_to_wv(theta_value, sig2, tau);
    }
    // MA1 
    else if(element_type == "MA1"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      wv_theo.col(i) = ma1_to_wv(theta_value, sig2, tau);
    }
    // WN
    else if(element_type == "WN"){
      wv_theo.col(i) = wn_to_wv(theta_value, tau);
    }// DR
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
    // ARMA11
    else if(element_type == "ARMA11"){
      ++i_theta;
      double th = theta(i_theta);
    
      ++i_theta;
      double sig2 = theta(i_theta);
      
      wv_theo.col(i) = arma11_to_wv(theta_value, th, sig2, tau);
    }
    else { // "SARIMA"
      
      // Unpackage ARMA model parameter
      arma::vec model_params = objdesc(i);
      
      unsigned int pop = arma::sum(model_params.rows(0,3));
      
      // Extract theta_values (includes the active theta_value)
      // Takes np + nq + nsp + nsq values
      arma::vec theta_values = theta.rows(i_theta, i_theta + pop - 1);
      
      // Increase the theta counter
      i_theta += pop;
      
      // Setup parameters
      arma::field<arma::vec> psetup = sarma_expand(theta_values, model_params);
      
      wv_theo.col(i) = arma_to_wv(psetup(0), psetup(1), theta(i_theta), tau);
    }
    
    ++i_theta;
  }

  return wv_theo;
}

//' Decomposed WV to Single WV
//' 
//' This function computes the combined processes to WV (haar) in a given model.
//' @param decomp A \code{mat} with scales as rows and processes as columns
//' @return A \code{vec} containing the wavelet variance of the process for the overall model
//' @template misc/haar_wv_formulae_link
//' @examples
//' model = AR1(.3,2) + DR(.001)
//' ntau = 8
//' tau = 2^(1:ntau)
//' wv.theo = decomp_theoretical_wv(model$theta, model$desc, model$objdesc, tau)
//' wv.total = decomp_to_theo_wv(wv.theo)
//' @keywords internal
// [[Rcpp::export]]
arma::vec decomp_to_theo_wv(const arma::mat& decomp){

  return arma::sum(decomp, 1);
}

/* ------------------------------ End Process to WV Functions --------------------------------- */
