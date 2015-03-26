#include <RcppArmadillo.h>
using namespace Rcpp;

#include "transform_data.h"

//' @title Pseudo Logit Inverse Function
//' @description This function computes the pseudo inverse of a logit transformation of the parameters in order to constrain them to a positive domain 
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @examples
//' x.sim = rnorm(100)
//' pseudo_logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit_inv(const arma::vec& x){
  return 2*exp(x)/(1 + exp(x)) -1;
}

//' @title Logit Inverse Function
//' @description This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @examples
//' x.sim = rnorm(100)
//' logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec logit_inv(const arma::vec& x){
  return 1/(1 + exp(-x));
}

//' @title Pseudo Logit Function
//' @description This function compute the link function to constrain parameters to a positive domain.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @examples
//' x.sim = runif(100)
//' pseudo_logit(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit(const arma::vec& x){
  arma::vec p = (x+1)/2;
  return log(p/(1 - p));
}

double pseudo_logit(double x){
  double p = (x+1)/2;
  return log(p/(1 - p));
}

//' @title Logit Function
//' @description This function computes the logit link function.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @examples
//' x.sim = runif(100)
//' logit(x.sim)
// [[Rcpp::export]]
arma::vec logit(const arma::vec& x){
  return log(x/(1 - x));
}

double logit(double x){
  return log(x/(1 - x));
}