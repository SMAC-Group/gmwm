#include <RcppArmadillo.h>
#include "analytical_matrix_derivatives.h"
using namespace Rcpp;

// Numerical Derivatives

inline double square(double x){
  return x*x;
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
//' deriv_AR1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_AR1(double phi, double sig2, arma::vec tau){
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
     D.col(0) = (8.0*sig2)/(pow(phi - 1.0, 4.0)*pow(phi + 1.0, 2.0)*tausq) // common term (8v^2)/[(p-1)^4(p+1)^2*tau^2]
                % ( -2.0 * (phi * ( phi * (tau - 6.0) - 4.0) - tau - 2.0) % phi_tau_1ov2
                      + (phi * ( phi * (tau - 3.0) - 2.0) - tau - 1.0) % phi_tau
                      - phi * (phi * (phi * tau + tau + 9) - tau + 6)
                      + tau - 3.0
                  );
                  
     // partial derivative with respect to sig2
     D.col(1) = 4.0*(phi*(-8.0*phi_tau_1ov2+2.0*phi_tau+phi*tau+6.0)-tau)/(pow(phi-1.0,3.0)*(phi + 1.0)*tausq);
     
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
//' deriv_2nd_AR1(.3, 1, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_AR1(double phi, double sig2, arma::vec tau){
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
     D.col(0) = (8.0*sig2)/(pow(phi - 1.0, 5.0)*phi*pow(phi + 1.0, 3.0)*tausq) // common term (8v^2)/[(p-1)^5*p*(p+1)^3*tau^2]
                % ( -1 * (phi * ( phi * (phi * (tau - 8.0) % (phi * (tau - 6.0) - 8.0) - 2.0 * (tau - 6.0) % tau +64.0 ) + 8.0*(tau+2.0) ) + tau % (tau + 2.0) ) % phi_tau_1ov2
                      + (phi * (phi * (phi * (tau - 4.0) % (phi * (tau - 3.0) - 4.0) - 2.0 * (tau - 3.0) % tau + 16.0 ) + 4 * (tau + 1) ) + tau % (tau + 1.0) ) % phi_tau
                      + 3.0 * phi * ( phi * (phi * (square(phi)*tau + 2.0 * phi * (tau + 6.0) + 16.0 ) - 2.0 * (tau - 8.0) ) - tau + 4 )
                  );
                  
     // partial derivative with respect to sig2
     D.col(1).fill(0);
     
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
//' deriv_DR(5.3, 2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_DR(double omega, arma::vec tau){
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
//' deriv_2nd_DR(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_2nd_DR(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
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
//' deriv_QN(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_QN(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = (2.0*arma::square(tau)+1.0)/(24.0*tau);
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
//' deriv_RW(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_RW(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = (2.0*arma::square(tau)+1.0)/(24.0*tau);
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
//' deriv_WN(2^(1:5))
// [[Rcpp::export]]
arma::mat deriv_WN(arma::vec tau){
     unsigned int ntau = tau.n_elem;
     arma::mat D(ntau, 1);
     D.col(0) = 1.0/tau;
     return D;
}
