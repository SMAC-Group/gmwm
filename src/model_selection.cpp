#include <RcppArmadillo.h>

#include "model_selection.h"

using namespace Rcpp;

//' @title D Matrix
//' @description 
//' @param At_j A \code{mat} that is of dimensions J x P containing the sederivative of A(theta_hat_j)/dTheta_hat_(j,tau)
//' @param omega A \code{mat} that is of dimension P x P used in obtaining the GMWM estimator.
//' @param diff A \code{vec} that is the difference of the WV empirical and WV theoretical
//' @return 
//' @details
//' hi
// [[Rcpp::export]]
arma::mat D_matrix(const arma::mat& At_j, const arma::mat& omega, const arma::vec& diff){
  unsigned int p = At_j.n_cols;
  unsigned int J = At_j.n_rows;
  arma::vec save = omega*diff; // (J x J) * (J x 1) = J x 1
  arma::mat D(p, p);

  arma::mat temp = arma::zeros<arma::mat>(p, J);
  for(unsigned int i = 0; i < p; i++){
    // force to column vector
    temp.row(i) = arma::trans(At_j.col(i));
    
    D.col(i) = temp*save; // (P x J) * (J x 1) = P x 1
    temp.row(i).fill(0);
  }
  
  return D;
}

//' @title Model Score
//' @description Calculates the modeling score of a GMWM
//' @param At
//' @param At_j 
//' @param omega
//' @param v_hat
//' @param diff A \code{vec} that is the difference of the WV empirical and WV theoretical
//' @param T An \code{unsigned int} that is awesome!
//' @return 
//' @details
//' hi
// [[Rcpp::export]]
arma::mat B_matrix(arma::mat A, arma::mat a_omega){
  unsigned int p = A.n_cols;
  arma::mat B(p,p);
  
  for(unsigned int i = 0; i < p; i++){
    B.col(i) = a_omega* A.col(i);
  }
  
  return B;
}

//' @title Model Score
//' @description Calculates the modeling score of a GMWM
//' @param A A \code{mat} that contains the first derivatives of the processes
//' @param At_j  A \code{mat} that contains the second derivative of each process
//' @param omega A \code{mat} that contains the omega used when calculating the GMWM
//' @param v_hat A \code{mat} that contains the covariance matrix
//' @param diff A \code{vec} that is the difference of the WV empirical and WV theoretical
//' @param T An \code{unsigned int} that is awesome!
//' @return 
//' @details
//' hi
// [[Rcpp::export]]
arma::vec model_score(arma::mat A, arma::mat At_j, arma::mat omega, arma::mat v_hat, arma::vec diff, unsigned int N){
  
  arma::mat At = arma::trans(A);
  arma::mat D = D_matrix(At_j, omega, diff);
  arma::mat B = B_matrix(A, At*omega);
  
  arma::mat d_b = D-B;
  arma::mat db_t = arma::trans(d_b);
  
  arma::mat dTheta = -arma::inv(db_t * d_b)*db_t*At*omega;

  arma::vec score_info(3);
  
  score_info(0) = arma::as_scalar(arma::trans(diff)*omega*diff + 2.0/double(N)*arma::trace(A * dTheta * omega * v_hat));
  score_info(1) = arma::as_scalar(arma::trans(diff)*omega*diff);
  score_info(2) = 2.0/double(N)*arma::trace(A * dTheta * omega * v_hat);
  
  // Return the score
  return score_info;
}