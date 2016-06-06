#include <RcppArmadillo.h>

#include "lm.h"

//' MLR in Armadillo
//' 
//' Perform Multiple Linear Regression using armadillo in C++
//' 
//' @param y A \code{vec} of length \eqn{N\times 1}{N x 1} containing the responses.
//' @param X A \code{mat} with dimensions \eqn{N \times p}{N x p}, which is the design matrix.
//' @return A \code{field<vec>} with:
//' \describe{
//'   \item{coef}{Coefficients}
//'   \item{resid}{Residuals}
//'   \item{sigma2}{Sigma^2}
//' }
//' @keywords internal
//' @examples
//' x = cbind(1,1:10)
//' y = cumsum(rep(.45,10))
//' 
//' lm_arma(y, x)[[1]]
//' 
//' coef(lm(y~x-1))
// [[Rcpp::export]]
arma::field<arma::vec> lm_arma(const arma::vec & y, const arma::mat & X) {
  
  int n = X.n_rows, p = X.n_cols;
  
  arma::vec coef = arma::solve(X,y);  
  
  arma::vec resid = y - X*coef; 
  
  arma::vec sig2 = arma::trans(resid)*resid/(n-p);
  
  arma::field<arma::vec> o(3);
  
  o(0) = coef;
  o(1) = resid;
  o(2) = sig2;
  
  return o;
}

//' Linear Regression with Drift
//' 
//' Perform a linear regression with drift.
//' 
//' @param x A \code{vec} of length \eqn{N\times 1}{N x 1} containing the responses.
//' @return A \code{field<vec>} with:
//' \describe{
//'   \item{coef}{Coefficients}
//'   \item{resid}{Residuals}
//'   \item{sigma2}{Sigma^2}
//' }
//' @keywords internal
//' @examples
//' x = 1:10
//' y = cumsum(rep(.7,10))
//' 
//' lm_dr(y)[[1]]
//' 
//' coef(lm(y~x-1))
//' 
// [[Rcpp::export]]
arma::field<arma::vec> lm_dr(const arma::vec & x) {
  
  int n =  x.n_elem;
  
  // No need to remove 1 from N. 
  arma::vec t = arma::linspace<arma::vec>(1,  n,  n);

  return lm_arma(x, t);
}