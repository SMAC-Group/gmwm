/* Copyright (C) 2014 - 2016 James Balamuta, Stephane Guerrier, Roberto Molinari
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

#include "additional_moments.h"

/* ------------------------ START ADDITIONAL MOMENTS ------------------------ */

//' @title Expected value DR
//' @description This function computes the expected value of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the expected value of the drift.
//' @keywords internal
//' @examples
//' e_drift(1,200)
// [[Rcpp::export]]
arma::vec e_drift(double omega, int n_ts){
  arma::vec out(1);
  out(0) = omega*(n_ts + 1.0)/2.0;
  return out;
}

//' @title Second moment DR
//' @description This function computes the second moment of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the second moment of the drift.
//' @keywords internal
//' @examples
//' m2_drift(1, 200)
// [[Rcpp::export]]
arma::vec m2_drift(double omega, int n_ts){
  arma::vec out(1);
  out(0)=(omega*omega)*(double(n_ts*n_ts)/3.0 + double(n_ts)/2.0 + 1.0/6.0);
  return out;
}

//' @title Variance DR
//' @description This function computes the variance of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the variance of the drift.
//' @keywords internal
//' @examples
//' var_drift(1, 200)
// [[Rcpp::export]]
arma::vec var_drift(double omega, int n_ts){
  // Compute m1
  arma::vec m1 = e_drift(omega, n_ts);
	
	// Compute m2
	arma::vec m2 = m2_drift(omega, n_ts);
	
	// Compute var
  return (m2 - m1*m1)*double(n_ts)/double(n_ts-1.0);
}

/* ------------------------ END ADDITIONAL MOMENTS ------------------------ */

