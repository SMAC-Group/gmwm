#include <RcppArmadillo.h>

#include "model_selection.h"
using namespace Rcpp;

// [[Rcpp:export]]
arma::mat D_matrix(arma::mat At_j, arma::mat omega, arma::mat diff){
  unsigned int p = At_j.n_cols;
  unsigned int J = At_j.n_rows;
  arma::vec save = omega*diff;
  arma::mat D(p, p);
  
  arma::mat temp = arma::zeros<arma::mat>(J,p);
  for(unsigned int i = 0; i < p; i++){
    temp.col(i) = At_j.col(i);
    D.col(i) = temp*save;
    temp.col(i).fill(0);
  }
  
  return D;
}

// [[Rcpp:export]]
arma::mat B_matrix(arma::mat A, arma::mat a_omega){
  unsigned int p = A.n_cols;
  arma::mat B(p,p);
  
  for(unsigned int i = 0; i < p; i++){
    B.col(i) = a_omega* A.col(i);
  }
  
  return B;
}

//' @title
//' @param At
//' @param At_j 
//' @param omega
//' @param v_hat
//' @param diff A \code{vec} that is the difference of the WV empirical and WV theoretical
//' @param T An \code{unsigned int} that is awesome!
//' @return 
//' @details
//' hi
// [[Rcpp:export]]
double model_score(arma::mat A, arma::mat At_j, arma::mat omega, arma::mat v_hat, arma::vec diff, unsigned int T){
    
  arma::mat D = D_matrix(At_j, omega, diff);
  arma::mat B = B_matrix(A, A*omega);
  
  arma::mat d_b = D-B;
  arma::mat db_t = arma::trans(d_b);
  
  arma::mat dTheta = -arma::inv(db_t * d_b)*db_t*A*omega;
  
  double score = arma::as_scalar(arma::trans(diff)*omega*diff + 2.0/double(T)*arma::trace(At_j * dTheta * omega * v_hat));
  return score;
}