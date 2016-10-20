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

// Gen ARMA
#include "gen_process.h"

// WV
#include "wave_variance.h"

// MODWT
#include "dwt.h"

// HAAR FILTER
#include "wv_filters.h"

//' @title Indirect Inference for ARMA
//' @description Option for indirect inference
//' @param ar A \code{vec} that contains the coefficients of the AR process.
//' @param ma A \code{vec} that contains the coefficients of the MA process.
//' @param sigma2 A \code{double} that indicates the sigma2 parameter of the ARMA process.
//' @param N A \code{int} that indicates how long the time series is.
//' @param robust A \code{bool} that indicates whether the estimation should be robust or not.
//' @param eff A \code{double} that specifies the amount of efficiency required by the robust estimator.
//' @param H A \code{int} that indicates how many iterations should take place.
//' @return A \code{vec} with the indirect inference results.
//' @keywords internal
// [[Rcpp::export]]
arma::vec idf_arma(const arma::vec& ar, const arma::vec& ma,
                   const double sigma2,
                   unsigned int N, bool robust, double eff, unsigned int H){
  
  // Access seed functions
  Rcpp::Environment baseEnv("package:base");
  Rcpp::Function setSeed = baseEnv["set.seed"];
  
  setSeed(1);
  
  unsigned int nb_level = floor(log2(N));
  
  arma::mat id(nb_level, H);
  for(unsigned int i=0; i<H; i++){
    
    // Generate x_t ~ F_theta
    arma::vec x = gen_arma(N, ar, ma, sigma2, 0);
    
    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level, "periodic");
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter(), "modwt");
    
    // Obtain WV
    arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
    // Store WV into matrix
    id.col(i) = wv_x;
  }
  
  return mean(id,1);
}


//' @title Indirect Inference for ARMA
//' @description Option for indirect inference
//' @param ar A \code{vec} that contains the coefficients of the AR process.
//' @param ma A \code{vec} that contains the coefficients of the MA process.
//' @param sigma2 A \code{double} that indicates the sigma2 parameter of the ARMA process.
//' @param N A \code{int} that indicates how long the time series is.
//' @param robust A \code{bool} that indicates whether the estimation should be robust or not.
//' @param eff A \code{double} that specifies the amount of efficiency required by the robust estimator.
//' @param H A \code{int} that indicates how many iterations should take place.
//' @return A \code{vec} with the indirect inference results.
//' @keywords internal
// [[Rcpp::export]]
arma::vec idf_arma_total(const arma::vec& ar, const arma::vec& ma,
                    const double sigma2,
                    unsigned int N, bool robust, double eff, unsigned int H){
  
  // Access seed functions
  Rcpp::Environment baseEnv("package:base");
  Rcpp::Function setSeed = baseEnv["set.seed"];
  
  setSeed(1);
  
  unsigned int n = N*H;
  
  unsigned int nb_level = floor(log2(n));
  
  // Generate x_t ~ F_theta
  arma::vec x = gen_arma(n, ar, ma, sigma2, 0);
  
  // MODWT transform
  arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level, "periodic");
  arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter(), "modwt");
  arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
  return wv_x;
}
