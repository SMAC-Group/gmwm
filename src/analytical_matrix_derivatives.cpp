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

//' ARMA Adapter to ARMA to WV Process function
//' 
//' Molds the data so that it works with the arma_to_wv function.
//' @param theta A \code{vec} that contains all the parameter estimates.
//' @param p     A \code{int} that indicates the number of AR coefficients
//' @param q     A \code{int} that indicates the number of MA coefficients.
//' @template misc/tau
//' @return A \code{vec} containing the ARMA to WV results
//' @details 
//' The data conversion is more or less a rearrangement of values without using the obj desc. 
//' @keywords internal
//' @backref src/analytical_matrix_derivatives.cpp
//' @backref src/analytical_matrix_derivatives.h
//' @template author/jjb
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
  
  return arma_to_wv(ar, ma, sigma2, tau);
}


//' Calculates the Jacobian for the ARMA process
//' 
//' Take the numerical derivative for the first derivative of an ARMA using the 2 point rule.
//' @inheritParams arma_adapter
//' @return A \code{mat} that returns the first numerical derivative of the ARMA process.
//' @keywords internal
//' @backref src/analytical_matrix_derivatives.cpp
//' @backref src/analytical_matrix_derivatives.h
//' @template author/jjb
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
//' 
//' Obtain the first derivative of the AR(1) process. 
//' @param phi    A \code{double} corresponding to the phi coefficient of an AR(1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an AR(1) process.
//' @template misc/tau
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\phi}{phi} 
//' and the second column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2}
//' @template deriv_wv/1st/deriv1_ar1
//' @template author/jjb
//' @examples
//' deriv_ar1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_ar1(double phi, double sigma2, const arma::vec& tau){
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
     D.col(0) = (2.0*sigma2)/(pow(phi - 1.0, 4.0)*pow(phi + 1.0, 2.0)*tausq) // common term (8v^2)/[(p-1)^4(p+1)^2*tau^2]
                % (-3.0 + tau - 2.0 * phi_tau_1ov2 % ( - 2.0 - tau + phi*(-4 +(-6 + tau)*phi))
                      + phi_tau % (-1.0 - tau + phi* (-2.0 + (-3.0 + tau)*phi))  
                      - phi * (6.0 - tau + phi * (9.0 + tau + tau * phi))
                  );
                  
     // partial derivative with respect to sigma2
     D.col(1) = (phi*(-8.0*phi_tau_1ov2+2.0*phi_tau+phi*tau+6.0)-tau)/(pow(phi-1.0,3.0)*(phi + 1.0)*tausq);
     
     return D;
}

//' Analytic second derivative matrix for AR(1) process
//' 
//' Calculates the second derivative for the AR(1) process and places it into a matrix form.
//' The matrix form in this case is for convenience of the calculation. 
//' @param phi    A \code{double} corresponding to the phi coefficient of an AR(1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an AR(1) process.
//' @template misc/tau
//' @return A \code{matrix} with the first column containing the
//'  second partial derivative with respect to \eqn{\phi}{phi} and
//'   the second column contains the second partial derivative with 
//'   respect to \eqn{\sigma ^2}{sigma^2}
//' @template author/jjb
//' @examples
//' deriv_2nd_ar1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_ar1(double phi, double sigma2, const arma::vec& tau){
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
     D.col(0) = (2.0*sigma2)/(pow(phi - 1.0, 5.0)*phi*pow(phi + 1.0, 3.0)*tausq) // common term (8v^2)/[(p-1)^5*p*(p+1)^3*tau^2]
                % ( -1 * (phi * ( phi * (phi * (tau - 8.0) % (phi * (tau - 6.0) - 8.0) - 2.0 * (tau - 6.0) % tau +64.0 ) + 8.0*(tau+2.0) ) + tau % (tau + 2.0) ) % phi_tau_1ov2
                      + (phi * (phi * (phi * (tau - 4.0) % (phi * (tau - 3.0) - 4.0) - 2.0 * (tau - 3.0) % tau + 16.0 ) + 4 * (tau + 1) ) + tau % (tau + 1.0) ) % phi_tau
                      + 3.0 * phi * ( phi * (phi * (square(phi)*tau + 2.0 * phi * (tau + 6.0) + 16.0 ) - 2.0 * (tau - 8.0) ) - tau + 4 )
                  );

     // partial derivative with respect to phi and sigma2
     D.col(1) = (2.0/(tausq * pow(phi - 1.0, 4.0)*pow(phi + 1.0, 2.0))) % 
                  (-3.0 + tau - 2.0*phi_tau_1ov2 % (-2.0 - tau + phi * (-4.0 + (-6.0 +tau)*phi)) +
                    phi_tau % (-1.0 - tau +phi*(-2.0 + (-3.0 + tau)*phi)) - phi * (6.0 - tau + phi*(9.0+tau+phi*tau)) );

     // partial derivative with respect to sigma2
     D.col(2).fill(0);
     
     return D;
}

//' Analytic D matrix for MA(1) process
//' 
//' Obtain the first derivative of the MA(1) process. 
//' @param theta  A \code{double} corresponding to the theta coefficient of an MA(1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an MA(1) process.
//' @template misc/tau
//' @return A \code{matrix} with the first column containing the partial derivative with respect to \eqn{\theta}{theta}
//'  and the second column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2}
//' @template deriv_wv/1st/deriv1_ma1
//' @template author/jjb
//' @examples
//' deriv_ma1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_ma1(double theta, double sigma2, const arma::vec& tau){
  unsigned int ntau = tau.n_elem;
  arma::mat D(ntau,2);
  
  // sigma ^2/ tau^2
  arma::vec shared = (1 / arma::square(tau));
  
  // partial derivative with respect to theta
  D.col(0) = shared % (2.0 * (theta + 1.0) * tau - 6.0) * sigma2;
    
  // partial derivative with respect to sigma2
  D.col(1) = shared % (square(theta + 1.0) * tau - 6.0 * theta);
  
  return D;
}



//' Analytic second derivative for MA(1) process
//' 
//' To ease a later calculation, we place the result into a matrix structure. 
//' @param theta  A \code{double} corresponding to the theta coefficient of an MA(1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an MA(1) process.
//' @template misc/tau
//' @return A \code{matrix} with the first column containing the second partial derivative with respect to \eqn{\theta}{theta},
//'  the second column contains the partial derivative with respect to \eqn{\theta}{theta} and \eqn{\sigma ^2}{sigma^2},
//'  and lastly we have the second partial derivative with respect to \eqn{\sigma ^2}{sigma^2}.
//' @template author/jjb
//' @examples
//' deriv_2nd_ma1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_ma1(double theta, double sigma2, const arma::vec& tau){
  unsigned int ntau = tau.n_elem;
  arma::mat D(ntau,3);
  
  // partial derivative with respect to phi
  D.col(0) = 2.0*sigma2 / tau;
  
  // partial derivative with respect to phi and sigma2
  D.col(1) = sigma2 / pow(tau,4.0) % // shared sigma2 / tau^4
              (2.0 * (theta + 1.0) * tau - 6.0) % // (2*(theta+1)*tau - 6)
              (square(theta + 1.0) * tau - 6.0 * theta);
  
  // partial derivative with respect to sigma2
  D.col(2).fill(0);
  
  return D;
}

//' Analytic D matrix for Drift (DR) Process
//' 
//' Obtain the first derivative of the Drift (DR) process. 
//' @param omega A \code{double} that is the slope of the drift.
//' @template misc/tau
//' @return A \code{matrix} with the first column containing the partial derivative 
//' with respect to \eqn{\omega}{omega}.
//' @template deriv_wv/1st/deriv1_dr
//' @template author/jjb
//' @examples
//' deriv_dr(5.3, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_dr(double omega, const arma::vec& tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau ,1);
     D.col(0) = omega*arma::square(tau) / 8.0;
     return D;
}

//' Analytic second derivative matrix for drift process
//' 
//' To ease a later calculation, we place the result into a matrix structure. 
//' @template misc/tau
//' @return A \code{matrix} with the first column containing 
//' the second partial derivative with respect to \eqn{\omega}{omega}.
//' @template author/jjb
//' @examples
//' deriv_2nd_dr(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_dr(const arma::vec& tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau,1);
     D.col(0) = arma::square(tau) / 8.0;
     return D;
}

//' Analytic D matrix for Quantization Noise (QN) Process
//' 
//' Obtain the first derivative of the Quantization Noise (QN) process. 
//' @template misc/tau
//' @return A \code{matrix} with the first column containing 
//' the partial derivative with respect to \eqn{Q^2}{Q^2}.
//' @template deriv_wv/1st/deriv1_qn
//' @template author/jjb
//' @examples
//' deriv_qn(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_qn(const arma::vec& tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = 6/arma::square(tau);
     return D;
}

//' Analytic D matrix Random Walk (RW) Process
//' 
//' Obtain the first derivative of the Random Walk (RW) process. 
//' @template misc/tau
//' @return A \code{matrix} with the first column containing
//'  the partial derivative with respect to \eqn{\gamma^2}{gamma^2}.
//' @template deriv_wv/1st/deriv1_rw
//' @template author/jjb
//' @examples
//' deriv_rw(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_rw(const arma::vec& tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = (arma::square(tau)+2.0)/(12.0*tau);
     return D;
}

//' Analytic D Matrix for a Gaussian White Noise (WN) Process
//' 
//' Obtain the first derivative of the Gaussian White Noise (WN) process. 
//' @template misc/tau
//' @return A \code{matrix} with the first column containing 
//' the partial derivative with respect to \eqn{\sigma^2}{sigma^2}.
//' @template deriv_wv/1st/deriv1_wn
//' @template author/jjb
//' @examples
//' deriv_wn(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_wn(const arma::vec& tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = 1.0/tau;
     return D;
}

//' Analytic D matrix of Processes
//' 
//' This function computes each process to WV (haar) in a given model.
//' @param theta   A \code{vec} containing the list of estimated parameters.
//' @param desc    A \code{vector<string>} containing a list of descriptors.
//' @param objdesc A \code{field<vec>} containing a list of object descriptors.
//' @template misc/tau
//' @return A \code{matrix} with the process derivatives going down the column
//' @details
//' Function returns the matrix effectively known as "D"
//' @template author/jjb
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
    if(element_type == "AR1" || element_type == "GM"){
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      D.cols(i_theta-1,i_theta) = deriv_ar1(theta_value, sig2, tau);
    } // MA 1
    else if(element_type == "MA1"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      D.cols(i_theta-1,i_theta) = deriv_ma1(theta_value, sig2, tau);
    }
    // WN
    else if(element_type == "WN"){
      D.col(i_theta) = deriv_wn(tau);
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
    // ARMA
    else {
      arma::vec o = objdesc(i);
      
      unsigned int p = o(0);
      unsigned int q = o(1);
      
      D.cols(i_theta, i_theta + p + q) = jacobian_arma(
        theta.rows(i_theta, i_theta + p + q),
        p, q, tau);
      
      i_theta += p + q;
    }
    
    ++i_theta;
  }

  return D;
}


//' Analytic D matrix of Processes
//' 
//' This function computes each process to WV (haar) in a given model.
//' @param theta     A \code{vec} containing the list of estimated parameters.
//' @param desc      A \code{vector<string>} containing a list of descriptors.
//' @param objdesc   A \code{field<vec>} containing a list of object descriptors.
//' @template misc/tau
//' @param omegadiff A \code{vec} that contains the result of Omega * (wv_empir - wv_theo)
//' @return A \code{matrix} with the process derivatives going down the column
//' @details
//' Function returns the matrix effectively known as "D"
//' @template author/jjb
//' @examples
//' # TBA
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
    if(element_type == "AR1" || element_type == "GM"){
      
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
      
    }else if(element_type == "MA1"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Return matrix with d/dtheta^2, d/dthetasig, d/dsig^2
      arma::mat s = deriv_2nd_ma1(theta_value, sig2, tau).t();
      
      // Modify p columns for theta
      A_i.rows(i_theta-1, i_theta) = s.rows(0,1);
      
      // Calculate the D matrix column value for the theta
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
