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


//' Analytic D matrix for ARMA(1,1) process
//' 
//' Obtain the first derivative of the ARMA(1,1) process. 
//' @param phi    A \code{double} corresponding to the phi coefficient of an ARMA(1,1) process.
//' @param theta  A \code{double} corresponding to the theta coefficient of an ARMA(1,1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an ARMA(1,1) process.
//' @template misc/tau
//' @return A \code{matrix} with:
//' \itemize{
//' \item The \strong{first} column containing the partial derivative with respect to \eqn{\phi}{phi};
//' \item The \strong{second} column containing the partial derivative with respect to \eqn{\theta}{theta};
//' \item The \strong{third} column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2}.
//' }
//' @template deriv_wv/1st/deriv1_arma11
//' @template author/jjb
//' @examples
//' deriv_arma11(.3, .4, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_arma11(double phi, double theta, double sigma2, const arma::vec& tau){
  unsigned int ntau = tau.n_elem;
  arma::mat D(ntau, 3);

  arma::vec phi_tau(ntau);
  arma::vec phi_tau_1ov2(ntau);
  arma::vec phi_tau_1ov2_m1(ntau);
  arma::vec phi_tau_m1(ntau);
  
  arma::vec tausq = arma::square(tau);
  
  for(unsigned int i = 0; i<ntau; i++){
    phi_tau(i) = pow(phi, tau(i)); // phi^(tau)
    phi_tau_1ov2(i) = pow(phi, tau(i)/2.0); // phi^(tau/2)
    phi_tau_1ov2_m1(i) = pow(phi, tau(i)/2.0 - 1.0); // phi^(tau/2 - 1)
    phi_tau_m1(i) = pow(phi, tau(i) - 1.0); // phi^(tau - 1)
  }
  
  // partial derivative with respect to phi
  D.col(0) = 2.0*sigma2*((-(3.0 - 4.0*phi_tau_1ov2 + phi_tau))*(1 + phi*(2 + 3*phi) + square(theta)*(1 + phi*(2 + 3*phi)) + 2*theta*(1 + phi*(3 + phi + square(phi)))) + 
    ((-square(1.0 + theta))*(-1.0 + phi)*square(1 + phi) - 
    2.0*phi_tau_1ov2_m1*(theta + phi)*
    (1.0 + theta*phi)*(-1.0 + square(phi)) + 
    phi_tau_m1*(theta + phi)*(1.0 + theta*phi)*(-1.0 + square(phi))) % tau)/
    (std::pow(-1.0 + phi,4.0)*square(1.0 + phi)*tausq);
  
  // partial derivative with respect to theta
  D.col(1) = (2.0*sigma2*((1.0 + 2.0*theta*phi + square(phi))*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) +
    (1.0 + theta)*(-1.0 + square(phi))*tau)) / 
    (std::pow(-1.0 + phi,3.0)*(1.0 + phi)*tausq);
  
  // partial derivative with respect to sigma2
  D.col(2) = ((-2*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi_tau_1ov2 + phi_tau) - 
                    0.5*square(1 + theta)*(-1 + square(phi))*tau))/
                (std::pow((-1.0 + phi),3.0)*(1 + phi)*tausq));
  
  return D;
}



//' Analytic D matrix for ARMA(1,1) process
//' 
//' Obtain the second derivative of the ARMA(1,1) process. 
//' @param phi    A \code{double} corresponding to the phi coefficient of an ARMA(1,1) process.
//' @param theta  A \code{double} corresponding to the theta coefficient of an ARMA(1,1) process.
//' @param sigma2 A \code{double} corresponding to the error term of an ARMA(1,1) process.
//' @template misc/tau
//' @return A \code{matrix} with:
//' \itemize{
//' \item The \strong{first} column containing the second partial derivative with respect to \eqn{\phi}{phi};
//' \item The \strong{second} column containing the second partial derivative with respect to \eqn{\theta}{theta};
//' \item The \strong{third} column contains the second partial derivative with respect to \eqn{\sigma ^2}{sigma^2}.
//' \item The \strong{fourth} column contains the partial derivative with respect to \eqn{\phi}{phi} and \eqn{\theta}{theta}.
//' \item The \strong{fiveth} column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2} and \eqn{\phi}{phi}.
//' \item The \strong{sixth} column contains the partial derivative with respect to \eqn{\sigma ^2}{sigma^2} and \eqn{\theta}{theta}.
//' }
//' @template deriv_wv/2nd/deriv2_arma11
//' @template author/jjb
//' @examples
//' deriv_2nd_arma11(.3, .4, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_arma11(double phi, double theta, double sigma2, const arma::vec& tau){
  unsigned int ntau = tau.n_elem;
  arma::mat D(ntau, 6);
  
  arma::vec phi_tau(ntau);
  arma::vec phi_tau_1ov2(ntau);
  arma::vec phi_tau_1ov2_m1(ntau);
  arma::vec phi_tau_1ov2_m2(ntau);
  arma::vec phi_tau_m1(ntau);
  
  arma::vec tausq = arma::square(tau);
  
  for(unsigned int i = 0; i<ntau; i++){
    phi_tau(i) = pow(phi, tau(i)); // phi^(tau)
    phi_tau_1ov2(i) = pow(phi, tau(i)/2.0); // phi^(tau/2)
    phi_tau_1ov2_m1(i) = pow(phi, tau(i)/2.0 - 1.0); // phi^(tau/2 - 1)
    phi_tau_1ov2_m2(i) = pow(phi, tau(i)/2.0 - 2.0); // phi^(tau/2 - 2)
    phi_tau_m1(i) = pow(phi, tau(i) - 1.0); // phi^(tau - 1)
  }
  
  // partial second derivative with respect to phi
  D.col(0) = 
    (2.0*sigma2*(-12.0*square(1.0 + phi)*(
        (-(theta + phi))*(1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 +  phi_tau) - (0.5)*square(1+theta)*(-1.0 + square(phi))*tau) + 
         square(-1.0 + phi)*(-2.0*square(-1.0+theta)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) 
         + phi_tau_1ov2_m2%(-2.0 + phi_tau_1ov2)%tau*(-1.0 + square(phi))*(-phi - square(theta)*phi + theta*(1.0 + 4.0*phi + square(phi)))
        + phi_tau_1ov2_m2%(-1.0 + phi_tau_1ov2)%tausq*square(1.0 + phi)*(theta + phi + square(theta)*phi + theta*square(phi))) + 6.0*(-1.0 + phi)*(1.0 + phi)*
    ((theta + phi)*(1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) + (0.5)*square(1.0+theta)*(-1.0 + square(phi))*tau + (1.0 + phi)*((-theta)*(theta + phi)*
    (3.0 - 4.0*phi_tau_1ov2 + phi_tau) - (1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) - square(1.0+theta)*phi*tau - phi_tau_1ov2_m1%(-2.0 + phi_tau_1ov2)%tau*(theta + phi)*(1.0 + theta*phi))))) / (std::pow(-1.0 + phi,5.0)*std::pow(1.0 + phi,3.0)*tausq);
  
  // partial second derivative with respect to theta
  D.col(1) = (2.0*sigma2*(2.0*phi*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) + (-1.0 + square(phi))*tau))/(std::pow(-1.0 + phi,3.0)*(1.0 + phi)*tausq);
  
  // partial second derivative with respect to sigma2
  D.col(2).fill(0);
  
  // partial second derivative with respect to theta and phi
  D.col(3) = -2.0*sigma2*(2.0*(3.0 - 4*phi_tau_1ov2 + phi_tau)*(1.0 + phi*(3.0 + phi + square(phi)) + theta*(1.0 + phi*(2.0 + 3.0*phi)))
                            + (2.0*(1.0 + theta)*(-1.0 + phi)*square(1.0 + phi) + 2.0*phi_tau_1ov2_m1*(-1.0 + square(phi))*(1.0 + 2.0*theta*phi + square(phi)) 
                                 - phi_tau_m1*(-1.0 + square(phi))*(1.0 + 2.0*theta*phi + square(phi)))%tau)/((std::pow(-1.0 + phi,4.0)*square(1.0 + phi)*tausq));


  // partial second derivative with respect to sigma2 and phi
  D.col(4) = (2.0*((-1.0 + phi)*((-(theta + phi))*(1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) - (0.5)*square(1.0 + theta)*(-1.0 + square(phi))*tau) +
    3.0*(1.0 + phi)*((-(theta + phi))*(1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) - (0.5)*square(1.0 + theta)*(-1.0 + square(phi))*tau) - 
    (-1.0 + phi)*(1.0 + phi)*((-theta)*(theta + phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) -
    (1.0 + theta*phi)*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) - square(1.0 + theta)*phi*tau - phi_tau_1ov2_m1%(-2.0 + phi_tau_1ov2)%tau*(theta + phi)*(1 + theta*phi))))/((std::pow(-1.0 + phi,4.0)*square(1.0 + phi)*tausq));
  
  // partial second derivative with respect to sigma2 and theta
  D.col(5) = (2.0*((theta + 1.0)*(square(phi) - 1.0)*tau + (2.0*theta*phi + square(phi) + 1.0)*(phi_tau - 4.0*phi_tau_1ov2 + 3.0)))/
  (std::pow(phi - 1.0,3.0)*(phi + 1.0)*tausq);
  
  
  return D;
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
//' @template deriv_wv/2nd/deriv2_ar1
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
     arma::vec phi_tau_mi(ntau);
     arma::vec phi_tau_1ov2_mi(ntau);
     
     arma::vec tausq = arma::square(tau);
     
     for(unsigned int i = 0; i<ntau; i++){
         phi_tau(i) = pow(phi, tau(i)); // phi^(tau)
         phi_tau_1ov2(i) = pow(phi, tau(i)/2.0); // phi^(tau/2)
         phi_tau_mi(i) = pow(phi, tau(i)-1.0); // phi^(tau-1)
         phi_tau_1ov2_mi(i) = pow(phi, tau(i)/2.0-1.0); // phi^(tau/2 -1)
     }

     // partial derivative with respect to phi
     D.col(0) = (2*sigma2*(4.0*(1.0 + 3.0*phi)*(1.0 + phi + square(phi))*(3.0 - 4.0*phi_tau_1ov2 + phi_tau) +
                          (-1.0 + square(phi))*(3.0*square(1.0 + phi) + 2.0*phi_tau_1ov2_mi*(1 + phi*(4 + 7*phi)) - phi_tau_mi*(1 + phi*(4 + 7*phi)))%tau 
                          + phi_tau_1ov2_mi % (-1.0 + phi_tau_1ov2) % tausq * square(-1.0 + square(phi)))) / (std::pow(-1.0 + phi,5.0)*std::pow(1.0 + phi,3.0)*tausq);

     // partial derivative with respect to phi and sigma2
     D.col(1) = (2.0*((-(3.0 - 4.0*phi_tau_1ov2 + phi_tau))*(1.0 + phi*(2.0 + 3.0*phi)) 
                      + (-1.0 + square(phi))*(-1 - phi - 2*phi_tau_1ov2 + phi_tau)%tau)) /
                    (std::pow(-1.0 + phi,4.0)*square(1.0 + phi)*tausq);
     
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
//' @template deriv_wv/2nd/deriv2_ma1
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
  D.col(1) = (-6.0 + 2.0*(1.0 + theta)*tau)/arma::square(tau);
  
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
    else if(element_type == "ARMA11"){
      ++i_theta;
      
      double th = theta(i_theta);
      
      ++i_theta;
      
      double sig2 = theta(i_theta);
      
      // Compute theoretical WV
      D.cols(i_theta-2,i_theta) = deriv_arma11(theta_value, th, sig2, tau);
      
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
  
  // Number of Models to iterate over
  unsigned int num_desc = desc.size();
  
  // Number of Paramters
  unsigned int p = theta.n_elem;
  
  // Number of Scales
  unsigned int ntau = tau.n_elem;
  
  // A_i matrix to calculate D
  // P x J
  arma::mat A_i = arma::zeros<arma::mat>(p, ntau); 
  
  // The D matrix, that we hope is accurate.
  // P x P
  arma::mat D = arma::zeros<arma::mat>(p, p);
  
  // Begin the process of iterating over the models. 
  unsigned int i_theta = 0;
  for(unsigned int i = 0; i < num_desc; i++){

    // Pop the first element. 
    
    double theta_value = theta(i_theta);
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1" || element_type == "GM"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Return matrix with d/dphi^2, d/dphisigma2, d/dsigma2^2
      arma::mat s = deriv_2nd_ar1(theta_value, sig2, tau).t();
      
      // Modify p columns for phi
      // We have dphi^2 vs. dphisigma2
      A_i.rows(i_theta-1, i_theta) = s.rows(0,1);
      
      // Calculate the D matrix column value for the phi
      D.col(i_theta-1) = A_i * omegadiff;
      
      // Modify p columns for sig
      // We have dphisigma2 vs. dsigma^4
      A_i.rows(i_theta-1, i_theta) = s.rows(1,2);
      
      // Calculate the D matrix column value for the sigma
      D.col(i_theta) = A_i * omegadiff;
      
      // Clear the Ai matrix
      // Only need to remove the dphisigma2 row as the other is just 0!
      A_i.row(i_theta-1).fill(0);
      
    }else if(element_type == "MA1"){
      
      ++i_theta;
      double sig2 = theta(i_theta);
      
      // Return matrix with d/dtheta^2, d/dthetasigma2, d/dsigma2^2
      arma::mat s = deriv_2nd_ma1(theta_value, sig2, tau).t();
      
      // Modify p columns for theta
      // We have dtheta^2 vs. dthetasigma2
      A_i.rows(i_theta-1, i_theta) = s.rows(0,1);
      
      // Calculate the D matrix column value for the theta
      D.col(i_theta-1) = A_i * omegadiff;
      
      // Modify p columns for sig
      // We have dthetasigma2 vs. dsigma^4
      A_i.rows(i_theta-1, i_theta) = s.rows(1,2);
      
      // Calculate the D matrix column value for the sigma
      D.col(i_theta) = A_i * omegadiff;
      
      // Clear the Ai matrix
      // Only need to remove the dthetasigma2 row as the other is just 0!
      A_i.row(i_theta-1).fill(0);
      
    }
    else if(element_type == "ARMA11"){
      
      // We receive the phi value in theta_value 
      // I know, confusing... Don't blame me, I wrote it in the past (er present?)
      
      unsigned int start_i_theta = i_theta;
      
      ++i_theta;
      // theta value
      double th = theta(i_theta);
      
      ++i_theta;
      // sigma2
      double sig2 = theta(i_theta);
      
      // Return matrix with d/dtheta^2, d/dthetasig, d/dsig^2
      arma::mat s = deriv_2nd_arma11(theta_value, th, sig2, tau).t();
      
      /* We receive from 2nd arma11 the following:
       * 
       * phi^2
       * theta^2
       * sigma^4
       * phitheta
       * sigphi
       * sigtheta
       */
      
      // To ease the sufferring, let's create a temp matrix
      arma::mat temp(3,ntau);
      
      // ---- Begin to fill for phi
      
      // Derivatives w.r.t to phi
      
      // Fill with phi^2
      temp.row(0) = s.row(0);
      // Fill with dphitheta
      temp.row(1) = s.row(3);
      // Fill with dphisigma2
      temp.row(2) = s.row(4);
      
      
      // Modify p columns for theta
      // We have dphi^2, dphitheta, dphisigma2
      A_i.rows(start_i_theta, i_theta) = temp;
      
      // Calculate the D matrix column value for the theta
      D.col(start_i_theta) = A_i * omegadiff;
      
      // ---- Begin to fill for theta
      
      // Derivatives w.r.t to theta
      
      // Fill with dphitheta
      temp.row(0) = s.row(3);
      // Fill with dtheta2
      temp.row(1) = s.row(1);
      // Fill with dthetasigma2
      temp.row(2) = s.row(5);
      
      // Modify p columns for sig
      // We have dthetasigma2 vs. dsigma^4
      A_i.rows(start_i_theta, i_theta) = temp;
      
      // Calculate the D matrix column value for the sigma
      D.col(start_i_theta+1) = A_i * omegadiff;
      
      // ---- Begin to fill for sigma2
      
      // Derivatives w.r.t to sigma2
      
      // Fill with dphisigma
      temp.row(0) = s.row(4);
      // Fill with dthetasigma
      temp.row(1) = s.row(5);
      // Fill with dsigma2^2
      temp.row(2) = s.row(2);
      
      // Modify p columns for sig
      // We have dthetasigma2 vs. dsigma^4
      A_i.rows(start_i_theta, i_theta) = temp;
      
      // Calculate the D matrix column value for the sigma
      D.col(start_i_theta+2) = A_i * omegadiff;
      
      
      // Clear the Ai matrix for next term
      // Only need to remove the dthetasigma2 row as the other is just 0!
      A_i.rows(start_i_theta, start_i_theta+1).fill(0);
      
      
    }
    // DR
    else if(element_type == "DR"){
      
      A_i.row(i_theta) = deriv_2nd_dr(tau).t(); 
      D.col(i_theta) = A_i * omegadiff;
      
      A_i.row(i_theta).fill(0);
    }
    else{
      // Already zero! (YAYAYA!)
      
      // Or generic SARIMA that is not supported.
    }
    
    ++i_theta;
  }
  
  return D;
}
