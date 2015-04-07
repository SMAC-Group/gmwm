/*#include <RcppArmadillo.h>

#include "analytical_matrix_derivatives.h"

//#include "inference.h"
using namespace Rcpp;

//' @title Calculate the Psi matrix
//' @description Computes the Psi matrix using supplied parameters
//' @param D
//' @param v_hat
//' @param omega
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat calculate_psi_matrix(const arma::mat& D, const arma::mat& v_hat, const arma::mat& omega){  
  arma::mat B = arma::inv(arma::trans(D)*omega*D)*D*omega;
  
  return B*v_hat*arma::trans(B);
}

//' @title Generate the Confidence Interval for Theta Estimates
//' @description yaya
//' @param theta
//' @param psi
//' @param p
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat theta_ci(arma::vec theta, arma::mat psi, double alpha = 0.05){
  unsigned int nparams = theta.n_elem;
  arma::mat ci(nparams, 3);
  double z = R::qnorm(1-alpha, 0, 1, true, false);
  
  arma::vec se = sqrt(diagvec(psi));
  
  ci.col(0) = theta - z*se;
  ci.col(1) = theta + z*se;
  ci.col(2) = se;

  return ci;
}

//' @title Compute the GOF Test
//' @description yaya
//' @param theta
//' @param psi
//' @param p
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::vec gof_test(arma::mat psi, arma::vec wv_empir, arma::vec wv_theo, unsigned int nparams, unsigned int N){
  
  
  return out;
}



// [[Rcpp::export]]
arma::field<arma::mat> inference_summary(arma::vec theta, arma::mat v_hat, arma::mat omega, arma::vec wv_empir, double alpha = 0.05){
  
  
  
  // calculate delta
  arma::vec store = gof_test(psi, wv_empir, wv_theo, theta.n_elem, N);
  arma::mat ci = theta_ci(theta, psi, alpha);

  arma::field<arma::mat> out(2);

  
  return out;
} */