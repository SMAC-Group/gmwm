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
#include "analytical_matrix_derivatives.h"
#include "inline_functions.h"
#include "process_to_wv.h"
#include "rtoarmadillo.h"

using namespace Rcpp;

//' @title ARMA Adapter to ARMA to WV Process function
//' @description Molds the data so that it works with the arma_to_wv function.
//' @param theta A \code{vec} that contains all the parameter estimates.
//' @param p A \code{int} that indicates the number of AR coefficients
//' @param q A \code{int} that indicates the number of MA coefficients.
//' @param tau A \code{vec} that lists the scales of the process e.g. 2^(1:J)
//' @return A \code{vec} containing the ARMA to WV results
//' @keywords internal
//' @backref src/analytical_matrix_derivatives.cpp
//' @backref src/analytical_matrix_derivatives.h
// [[Rcpp::export]]
arma::vec arma_adapter(const arma::vec& theta,
                       unsigned int p,
                       unsigned int q,
                       const arma::vec& tau){
  
  arma::vec ar(p), ma(q);
  
  double sigma2;
  
  if(p != 0){
    ar = theta.rows(0, p - 1);
  }else{
    ar = arma::zeros<arma::vec>(p);
  }
  
  if(q != 0 ){
    ma = theta.rows(p, p+q-1);
  }else{
    ma = arma::zeros<arma::vec>(q);
  }
  
  sigma2 = theta(p+q);
  
  return arma_to_wv(ar, ma, tau, sigma2);
}


//' @title Calculates the Jacobian for the ARMA process
//' @description Figure out the Jacobian for an ARMA process.
//' @param theta A \code{vec} that contains all the parameter estimates.
//' @param p A \code{int} that indicates the number of AR coefficients
//' @param q A \code{int} that indicates the number of MA coefficients.
//' @param tau A \code{vec} that lists the scales of the process e.g. 2^(1:J)
//' @return A \code{mat} that returns the numerical jacobian of the ARMA process.
//' @keywords internal
//' @backref src/analytical_matrix_derivatives.cpp
//' @backref src/analytical_matrix_derivatives.h
// [[Rcpp::export]]
arma::mat jacobian_arma(const arma::vec& theta,
                        unsigned int p,
                        unsigned int q,
                        const arma::vec& tau){
  
  
  unsigned int n = theta.n_elem;
  
  unsigned int out_length = tau.n_elem;
  
  // Args
  double eps = .0001, d = .0001, zero_tol = sqrt(DBL_EPSILON/.0000007), // 7e-7
    r = 4, v = 2;
  
  
  arma::vec h = abs(d * theta) + eps * (abs(theta) < zero_tol);

  arma::mat out(out_length, n);
  
  // There has to be a better way to initialize a field of matrices...
  arma::field<arma::mat> st(n);
  arma::mat submat(out_length, r);

  for(unsigned int i = 0; i < n; i++){
    st(i) = submat;
  }
  
  arma::vec seq_n = seq_len_cpp(n) - 1;
  
  for (unsigned int k = 0; k < r; k++) {
    for (unsigned int i = 0; i < n; i++) {
      
      // We need to replace this badly with something more generic.
      st(i).col(k) = (arma_adapter(theta + h % (i == seq_n), p,q, tau) 
                        - arma_adapter(theta - h % (i == seq_n), p,q, tau)) / (2 * h(i));
      
    }
    h = h/v;
  }
  
  for (unsigned int i = 0; i < n; i++){
    for (unsigned int m = 1; m <= r - 1; m++){
      arma::mat act = st(i);

      st(i) = ( act.cols(1, r - m) * pow(4,m) - act.cols(0,(r - m - 1)) )/(pow(4,m) - 1);
    }
    
    out.col(i) = st(i);
  }
  
  return out;
}


//' Analytic D matrix for AR(1) process
//' @param phi A \code{double} corresponding to the phi coefficient of an AR(1) process.
//' @param sig2 A \code{double} corresponding to the error term of an AR(1) process.
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\phi ^2}{sigma^2} and the second column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2}
//' @details
//' The haar wavelet variance is given as \eqn{\frac{{\left( {\frac{\tau }{2} - 3{\rho _0} - \frac{{\tau \rho _0^2}}{2} + 4\rho _0^{\frac{\tau }{2} + 1} - \rho _0^{\tau  + 1}} \right)\nu _0^2}}{{\frac{{{\tau ^2}}}{8}{{\left( {1 - {\rho _0}} \right)}^2}\left( {1 - \rho _0^2} \right)}}}{See PDF Manual for equation}
//' Note: \eqn{\phi = \rho}{phi = rho} and \eqn{V _0^2 = \sigma _0^2}{V[0]^2 = sigma[0]^2}.
//' Due to length, the analytical derivations of the AR(1) haar wavelet variance are given in a supplied file within vignette.
//' @author JJB
//' @examples
//' deriv_ar1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_ar1(double phi, double sig2, arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau,2);
     
     //double phi_n = pow(phi, 0.5); // phi^(1/2)
     
     arma::vec phi_tau(ntau);
     arma::vec phi_tau_1ov2(ntau);
     arma::vec tausq = arma::square(tau);
     
     for(unsigned int i = 0; i<ntau; i++){
         phi_tau(i) = pow(phi, tau(i)); // phi^(tau)
         phi_tau_1ov2(i) = pow(phi, tau(i)/2.0); // phi^(tau/2)
     }

     // partial derivative with respect to phi
     D.col(0) = (2.0*sig2)/(pow(phi - 1.0, 4.0)*pow(phi + 1.0, 2.0)*tausq) // common term (8v^2)/[(p-1)^4(p+1)^2*tau^2]
                % (-3.0 + tau - 2.0 * phi_tau_1ov2 % ( - 2.0 - tau + phi*(-4 +(-6 + tau)*phi))
                      + phi_tau % (-1.0 - tau + phi* (-2.0 + (-3.0 + tau)*phi))  
                      - phi * (6.0 - tau + phi * (9.0 + tau + tau * phi))
                  );
                  
     // partial derivative with respect to sig2
     D.col(1) = (phi*(-8.0*phi_tau_1ov2+2.0*phi_tau+phi*tau+6.0)-tau)/(pow(phi-1.0,3.0)*(phi + 1.0)*tausq);
     
     return D;
}

//' Analytic second derivative matrix for AR(1) process
//' @param phi A \code{double} corresponding to the phi coefficient of an AR(1) process.
//' @param sig2 A \code{double} corresponding to the error term of an AR(1) process.
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the second partial derivative with respect to \eqn{\phi ^2}{sigma^2} and the second column contains the second partial derivative with respect to \eqn{\sigma ^2}{sigma^2}
//' @details
//' The haar wavelet variance is given as \eqn{\frac{{\left( {\frac{\tau }{2} - 3{\rho _0} - \frac{{\tau \rho _0^2}}{2} + 4\rho _0^{\frac{\tau }{2} + 1} - \rho _0^{\tau  + 1}} \right)\nu _0^2}}{{\frac{{{\tau ^2}}}{8}{{\left( {1 - {\rho _0}} \right)}^2}\left( {1 - \rho _0^2} \right)}}}{See PDF Manual for equation}
//' Note: \eqn{\phi = \rho}{phi = rho} and \eqn{V _0^2 = \sigma _0^2}{V[0]^2 = sigma[0]^2}.
//' Due to length, the analytical derivations of the AR(1) haar wavelet variance are given in a supplied file within vignette.
//' @author JJB
//' @examples
//' deriv_2nd_ar1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_ar1(double phi, double sig2, arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau,3);
     
     //double phi_n = pow(phi, 0.5); // phi^(1/2)
     
     arma::vec phi_tau(ntau);
     arma::vec phi_tau_1ov2(ntau);
     arma::vec tausq = arma::square(tau);
     
     for(unsigned int i = 0; i<ntau; i++){
         phi_tau(i) = pow(phi, tau(i)); // phi^(tau)
         phi_tau_1ov2(i) = pow(phi, tau(i)/2.0); // phi^(tau/2)
     }

     // partial derivative with respect to phi
     D.col(0) = (2.0*sig2)/(pow(phi - 1.0, 5.0)*phi*pow(phi + 1.0, 3.0)*tausq) // common term (8v^2)/[(p-1)^5*p*(p+1)^3*tau^2]
                % ( -1 * (phi * ( phi * (phi * (tau - 8.0) % (phi * (tau - 6.0) - 8.0) - 2.0 * (tau - 6.0) % tau +64.0 ) + 8.0*(tau+2.0) ) + tau % (tau + 2.0) ) % phi_tau_1ov2
                      + (phi * (phi * (phi * (tau - 4.0) % (phi * (tau - 3.0) - 4.0) - 2.0 * (tau - 3.0) % tau + 16.0 ) + 4 * (tau + 1) ) + tau % (tau + 1.0) ) % phi_tau
                      + 3.0 * phi * ( phi * (phi * (square(phi)*tau + 2.0 * phi * (tau + 6.0) + 16.0 ) - 2.0 * (tau - 8.0) ) - tau + 4 )
                  );

     // partial derivative with respect to phi and sig2
     D.col(1) = (2.0/(tausq * pow(phi - 1.0, 4.0)*pow(phi + 1.0, 2.0))) % 
                  (-3.0 + tau - 2.0*phi_tau_1ov2 % (-2.0 - tau + phi * (-4.0 + (-6.0 +tau)*phi)) +
                    phi_tau % (-1.0 - tau +phi*(-2.0 + (-3.0 + tau)*phi)) - phi * (6.0 - tau + phi*(9.0+tau+phi*tau)) );

     // partial derivative with respect to sig2
     D.col(2).fill(0);
     
     return D;
}


//' Analytic D matrix for drift process
//' @param omega A \code{double} that is the slope of the drift.
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\omega _0}{omega[0]}.
//' @details
//' The haar wavelet variance is given as \eqn{{\nu ^2}\left( \tau  \right) = \frac{{{\tau ^2}\omega _0^2}}{2}}{nu^2(tau) = tau^2 omega_0^2 / 2}.
//' Taking the derivative with respect to \eqn{\omega _0^2}{omega_0^2} yields: \eqn{\frac{\partial }{{\partial {\omega _0}}}{\nu ^2}\left( \tau  \right) = {\tau ^2}{\omega _0}}{tau^2 * omega_0}
//' @author JJB
//' @examples
//' deriv_dr(5.3, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_dr(double omega, arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau ,1);
     D.col(0) = omega*arma::square(tau);
     return D;
}

//' Analytic second derivative matrix for drift process
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the second partial derivative with respect to \eqn{\omega _0}{omega[0]}.
//' @details
//' The haar wavelet variance is given as \eqn{{\nu ^2}\left( \tau  \right) = \frac{{{\tau ^2}\omega _0^2}}{2}}{nu^2(tau) = tau^2 omega_0^2 / 2}.
//' Taking the derivative with respect to \eqn{\omega _0^2}{omega_0^2} yields: \eqn{\frac{\partial }{{\partial {\omega _0}}}{\nu ^2}\left( \tau  \right) = {\tau ^2}{\omega _0}}{tau^2 * omega_0}
//' Taking second derivative with respect to \eqn{\omega _0^2}{omega_0^2} yields: \eqn{\frac{{{\partial ^2}}}{{\partial \omega _0^2}}{\nu ^2}\left( \tau  \right) = {\tau ^2}}{tau^2}
//' @author JJB
//' @examples
//' deriv_2nd_dr(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_dr(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau,1);
     D.col(0) = arma::square(tau);
     return D;
}

//' Analytic D matrix quantisation noise process
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{Q _0^2}{Q[0]^2}.
//' @details
//' The haar wavelet variance is given as \eqn{{\nu ^2}\left( \tau  \right) = \frac{{3Q_0^2}}{{2{\tau ^2}}}}{nu^2(tau) = 3*Q[0]^2 / 2*tau^2}.
//' Taking the derivative with respect to \eqn{Q _0^2}{Q[0]^2} yields: \deqn{\frac{\partial }{{\partial Q_0^2}}{\nu ^2}\left( \tau  \right) = \frac{3}{{2{\tau ^2}}}}{3/(2*tau^2)}.
//' The second derivative derivative with respect to \eqn{Q _0^2}{Q[0]^2} is then: \deqn{\frac{{{\partial ^2}}}{{\partial Q_0^4}}{\nu ^2}\left( \tau  \right) = 0}{0}.
//' @author JJB
//' @examples
//' deriv_qn(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_qn(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = 3.0/(2.0*arma::square(tau));
     return D;
}

//' Analytic D matrix random walk process
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\gamma _0^2}{gamma[0]^2}.
//' @details
//' The haar wavelet variance is given as \eqn{{\nu ^2}\left( \tau  \right) = \frac{{\left( {2{\tau ^2} + 1} \right)\gamma _0^2}}{{24\tau }}}{nu^2(tau) = (2*tau^2+1)*gamma^2 / (24*tau)}.
//' Taking the first derivative with respect to \eqn{\gamma _0^2}{gamma_0^2} yields: \deqn{\frac{{{\partial ^2}}}{{\partial \gamma _0^4}}{\nu ^2}\left( \tau  \right) = 0}{(2*tau^2+1) / (24*tau)}
//' The second derivative derivative with respect to \eqn{\gamma _0^2}{gamma[0]^2} is then: \deqn{\frac{{{\partial ^2}}}{{\partial \sigma_0^4}}{\nu ^2}\left( \tau  \right) = 0}{0}.
//' @author JJB
//' @examples
//' deriv_rw(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_rw(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = (arma::square(tau)+2.0)/(12.0*tau);
     return D;
}

//' Analytic D matrix white noise process
//' @param tau A \code{vec} that contains the scales to be processed (e.g. 2^(1:J))
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\sigma _0^2}{sigma[0]^2}.
//' @details
//' The haar wavelet variance is given as \eqn{{\nu ^2}\left( \tau  \right) = \frac{{\sigma _0^2}}{\tau }}{nu^2(tau) = sigma_0^2 / tau}.
//' Taking the derivative with respect to \eqn{\sigma _0^2}{sigma_0^2} yields: \eqn{\frac{\partial }{{\partial \sigma _0^2}}{\nu ^2}\left( \tau  \right) = \frac{1}{\tau }}{1/tau}
//' @author JJB
//' @examples
//' deriv_wn(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_wn(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = 1.0/tau;
     return D;
}

//' Analytic D matrix of Processes
//' @description This function computes each process to WV (haar) in a given model.
//' @param theta A \code{vec} containing the list of estimated parameters.
//' @param desc A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @param tau A \code{vec} containing the scales e.g. 2^(1:J)
//' @return A \code{matrix} with the process derivatives going down the column
//' @details
//' Function returns the matrix effectively known as "D"
//' @author JJB
//' @examples
//' mod = AR1(.4,1) + WN(.2) + DR(.005)
//' derivative_first_matrix(mod$theta, mod$desc, mod$obj.desc, 2^(1:9))
// [[Rcpp::export]]
arma::mat derivative_first_matrix(const arma::vec& theta, 
                                  const std::vector<std::string>& desc,
                                  const arma::field<arma::vec>& objdesc,
                                  const arma::vec& tau){
                                  
  unsigned int num_desc = desc.size();
  arma::mat D = arma::zeros<arma::mat>(tau.n_elem, theta.n_elem);
    
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
      D.cols(i_theta-1,i_theta) = deriv_ar1(theta_value, sig2, tau);
    }
    else if(element_type == "ARMA"){
      arma::vec o = objdesc(i);
      
      unsigned int p = o(0);
      unsigned int q = o(1);
      
      D.cols(i_theta, i_theta + p + q) = jacobian_arma(
          theta.rows(i_theta, i_theta + p + q),
          p, q, tau);
      
      i_theta += p + q;
    }
    // DR
    else if(element_type == "DR"){
      D.col(i_theta) = deriv_dr(theta_value, tau);
    }
    // QN
    else if(element_type == "QN"){
      D.col(i_theta) = deriv_qn(tau);
    }
    // RW
    else if(element_type == "RW"){
      D.col(i_theta) = deriv_rw(tau);
    }
    // WN
    else{
      D.col(i_theta) = deriv_wn(tau);
    }
    
    ++i_theta;
  }

  return D;
}


//' Analytic D matrix of Processes
//' @description This function computes each process to WV (haar) in a given model.
//' @param theta A \code{vec} containing the list of estimated parameters.
//' @param desc A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @param tau A \code{vec} containing the scales e.g. 2^(1:J)
//' @param omegadiff A \code{vec} that contains the result of Omega * (wv_empir - wv_theo)
//' @return A \code{matrix} with the process derivatives going down the column
//' @details
//' Function returns the matrix effectively known as "D"
//' @author JJB
//' @examples
//' #TBA
//' @keywords internal
// [[Rcpp::export]]
arma::mat D_matrix(const arma::vec& theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   const arma::vec& tau, const arma::vec& omegadiff){
  
  unsigned int num_desc = desc.size();
  
  
  unsigned int p = theta.n_elem;
  
  unsigned int ntau = tau.n_elem;
  
  // P x J
  arma::mat A_i = arma::zeros<arma::mat>(p, ntau); 
  
  // P x P
  arma::mat D = arma::zeros<arma::mat>(p, p);
  
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){
    // Add ARMA
    
    double theta_value = theta(i_theta);
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Return matrix with d/dtheta^2, d/dthetasig, d/dsig^2
      arma::mat s = deriv_2nd_ar1(theta_value, sig2, tau).t();
      
      // Modify p columns for phi
      A_i.rows(i_theta-1, i_theta) = s.rows(0,1);
      
      // Calculate the D matrix column value for the phi
      D.col(i_theta-1) = A_i * omegadiff;
      
      // Modify p columns for sig
      A_i.rows(i_theta-1, i_theta) = s.rows(1,2);
      
      // Calculate the D matrix column value for the sigma
      D.col(i_theta) = A_i * omegadiff;
      
      // Clear the Ai matrix
      A_i.row(i_theta-1).fill(0);
      
    }
    else if(element_type == "ARMA"){
      // implement later      
      
    }
    // DR
    else if(element_type == "DR"){
      
      A_i.row(i_theta) = deriv_2nd_dr(tau).t(); 
      D.col(i_theta) = A_i * omegadiff;
      
      A_i.row(i_theta).fill(0);
    }
    else{
      // Already zero! (YAYAYA!)
    }
    
    ++i_theta;
  }
  
  return D;
}
