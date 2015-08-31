#include <RcppArmadillo.h>

#include "model_selection.h"

// Needed for COV bootstrap
#include "bootstrappers.h"

using namespace Rcpp;


/* To whom it may concern....
 * 
 * The D matrix is located in the analytical_matrix_derivatives.cpp code file.
 * 
 */

//' @title B Matrix
//' @description B Matrix
//' @param A A \code{mat} containing the first derivatives of the process.
//' @param at_omega A \code{mat} containing A^T * Omega
//' @return A \code{mat}
//' @details
//' TBA
// [[Rcpp::export]]
arma::mat B_matrix(const arma::mat& A, const arma::mat& at_omega){
  unsigned int p = A.n_cols;
  arma::mat B(p,p);
  
  for(unsigned int i = 0; i < p; i++){
    B.col(i) = at_omega* A.col(i);
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
//' @return A \code{vec}
//' @details
//' The equation is slightly different than that stated in the paper due to the bootstrap already incorporating in 
//' N.
// [[Rcpp::export]]
arma::vec model_score(arma::mat A, arma::mat D, arma::mat omega, arma::mat v_hat, double obj_value){
  
  arma::mat At = arma::trans(A);
  arma::mat B = B_matrix(A, At*omega);
  
  arma::mat d_b = D-B;
  arma::mat db_t = arma::trans(d_b);
  arma::mat dTheta = -arma::pinv(db_t * d_b)*db_t*At*omega;

  arma::vec score_info(2);
  
  score_info(0) = obj_value + arma::as_scalar(2.0*arma::trace(A * dTheta * omega * v_hat));
  score_info(1) = 2.0*arma::trace(A * dTheta * omega * v_hat);
  
  // Return the score
  return score_info;
}
