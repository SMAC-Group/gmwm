/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
 * (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#include <RcppArmadillo.h>

#include "arima_gmwm.h"

//' @title Hook into R's ARIMA function
//' @description Uses R's ARIMA function to obtain CSS values for starting condition
//' @param data A \code{vec} of data.
//' @param params A \code{vec} of the ARMA parameters
//' @return A \code{vec} containing the CSS of the ARMA parameters.
//' @keywords internal
// [[Rcpp::export]]
arma::vec Rcpp_ARIMA(const arma::vec& data,
                     const arma::vec& params){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function arima = stats["arima"];    
  
  arma::vec aparams(3);
  
  aparams(0) = params(0);
  aparams(1) = 0;
  aparams(2) = params(1);
  
  Rcpp::List Opt= arima(Rcpp::_["x"] = data,
                        Rcpp::_["order"] = aparams,
                        Rcpp::_["include.mean"] = false,
                        Rcpp::_["method"] = "CSS");
  
  arma::vec out = arma::join_cols(Rcpp::as<arma::vec>(Opt[0]), Rcpp::as<arma::vec>(Opt[1]));
  
  return out;
}
