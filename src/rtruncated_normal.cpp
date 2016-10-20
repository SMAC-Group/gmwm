/* Copyright (C) 2016  James Balamuta
 *
 * This file is part of GMWM R Methods Package
 *
 * The file uses methods in the rgen project and is free software: you can redistribute it and/or modify it
 * under the terms of the MIT License.
 *
 * The rgen project is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

#include <RcppArmadillo.h>

#include "rtruncated_normal.h"

//' Truncated Normal Distribution Sampling Algorithm
//' 
//' Enables sampling from a truncated normal
//' @param n      An \code{unsigned int} indicating the number of observations to generate.
//' @param mu     A \code{double} indicating the mean of the normal.
//' @param sigma  A \code{double} indicating the standard deviation of the normal.
//' @param a      A \code{double} that is the lower bound of the truncated normal.
//' @param b      A \code{double} that is the upper bound of the truncated normal.
//' @examples
//' rtruncated_normal(10, 1, 1, -2, 2)
// [[Rcpp::export]]
arma::vec rtruncated_normal(unsigned int n, double mu, double sigma, double a, double b){
  double phi_a = ( a - mu ) / sigma;
  double phi_b = ( b - mu ) / sigma;
  
  double phi_a_cdf = R::pnorm( phi_a, 0.0, 1.0, true, false );
  double phi_b_cdf = R::pnorm( phi_b, 0.0, 1.0, true, false ); 
  
  arma::vec o(n);
  for(unsigned int i = 0; i < n; i++){
    o(i) = sim_truncated_normal(phi_a_cdf, phi_b_cdf, mu, sigma);
  }
  
  // Release output
  return o;
}

double sim_truncated_normal(double phi_a_cdf, double phi_b_cdf, double mu, double sigma){
  // Generate u ~ U[0,1]
  double u = R::runif(0.0,1.0);
  
  double tnorm_value = phi_a_cdf + u * ( phi_b_cdf - phi_a_cdf );
  
  double tnorm_quantile = R::qnorm(tnorm_value, 0.0, 1.0, true, false);
  
  return mu + sigma * tnorm_quantile;
}

double rtruncated_normal(double mu, double sigma, double a, double b)
{
  double phi_a = ( a - mu ) / sigma;
  double phi_b = ( b - mu ) / sigma;
  
  double phi_a_cdf = R::pnorm( phi_a, 0.0, 1.0, true, false );
  double phi_b_cdf = R::pnorm( phi_b, 0.0, 1.0, true, false ); 
  
  return sim_truncated_normal(phi_a_cdf, phi_b_cdf, mu, sigma);
}
