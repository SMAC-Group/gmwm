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

#include "complex_tools.h"

//' @title Absolute Value or Modulus of a Complex Number Squared.
//' @description Computes the squared value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus squared for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod\left(z\right) = \sqrt{\left(x^2 + y^2\right)}}{r = Mod(z) = sqrt(x^2 + y^2)}
//' This function will return: \eqn{r^2 = Mod\left(z\right)^2 = x^2 + y^2}
//' @examples
//' Mod_squared_cpp(c(1+.5i, 2+1i, 5+9i))
//' @keywords internal
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
//' @keywords internal
// [[Rcpp::export]]
arma::vec Mod_cpp(const arma::cx_vec& x){
  return arma::sqrt(arma::square(arma::real(x)) + arma::square(arma::imag(x)));
}
