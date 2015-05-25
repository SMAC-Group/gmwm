#include <RcppArmadillo.h>

#include "complex_tools.h"
using namespace Rcpp;

//' @title Absolute Value or Modulus of a Complex Number Squared.
//' @description Computes the squared value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus squared for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod\left(z\right) = \sqrt{\left(x^2 + y^2\right)}}{r = Mod(z) = sqrt(x^2 + y^2)}
//' This function will return: \eqn{r^2 = Mod\left(z\right)^2 = x^2 + y^2}
//' @examples
//' Mod_squared_cpp(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_squared_cpp(const arma::cx_vec& x){
  return arma::square(arma::real(x)) + arma::square(arma::imag(x));
}

//' @title Absolute Value or Modulus of a Complex Number.
//' @description Computes the value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod(z) = \sqrt{(x^2 + y^2)}}
//' @examples
//' Mod_cpp(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_cpp(const arma::cx_vec& x){
  return arma::sqrt(arma::square(arma::real(x)) + arma::square(arma::imag(x)));
}
