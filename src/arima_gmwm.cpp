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
  
  // Bridge to R
  Rcpp::NumericVector aparams(3);
  
  aparams(0) = params(0);
  aparams(2) = params(1);
  
  // Load difference of 0.  (Taken care of in the gmwm_master)
  aparams(1) = 0; //params(6);
  
  Rcpp::List Opt;
  
  // Determine whether it is ARMA11 or SARIMA object
  // And check to see if there is seasonality (s)
  if(params.n_elem > 3 && params(5) > 0){
    
    Rcpp::NumericVector seasonal_order(3);
    
    // Seasonal AR
    seasonal_order(0) = params(2);
    
    // Seasonal Difference (Taken care of in the gmwm_master)
    seasonal_order(1) = 0; //params(7);
    
    // Seasonal MA
    seasonal_order(2) = params(3);
    
    Rcpp::List seasons = Rcpp::List::create(Rcpp::Named("order") = seasonal_order,
                                            Rcpp::Named("period") = params(5));

    Opt = arima(Rcpp::_["x"] = data,
                Rcpp::_["order"] = aparams,
                Rcpp::_["include.mean"] = false,
                Rcpp::_["seasonal"] = seasons,
                Rcpp::_["method"] = "CSS");
  }else{
    
    Opt = arima(Rcpp::_["x"] = data,
                Rcpp::_["order"] = aparams,
                Rcpp::_["include.mean"] = false,
                Rcpp::_["method"] = "CSS");
  }
  
  arma::vec out = arma::join_cols(Rcpp::as<arma::vec>(Opt[0]), Rcpp::as<arma::vec>(Opt[1]));
  
  return out;
}
