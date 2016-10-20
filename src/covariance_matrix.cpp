/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
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
#include "covariance_matrix.h"
#include "rtoarmadillo.h"

//' @title Computes the (MODWT) wavelet covariance matrix
//' @description Calculates the (MODWT) wavelet covariance matrix
//' @param signal_modwt A \code{field<vec>} that contains the modwt decomposition.
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @param compute_v A \code{string} that indicates what kind of matrix should be created. Possible options: "diag" or "none"
//' @param robust A \code{boolean} that triggers the use of the robust estimate.
//' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
//' @return A \code{field<mat>} containing the covariance matrix.
//' @keywords internal
//' @examples
//' \dontrun{
//' x=rnorm(100)
//' decomp = modwt(x)
//' V = compute_cov_cpp(decomp$data, decomp$nlevels, compute_v="diag", robust = TRUE, eff=0.6)
//' }
// [[Rcpp::export]]
arma::field<arma::mat> compute_cov_cpp(arma::field<arma::vec> signal_modwt, unsigned int nb_level, 
                                        std::string compute_v="diag",
                                        bool robust = true, double eff=0.6){// Compute asymptotic covariance matrix
  arma::field<arma::mat> V(2);
  //arma::mat V = diagmat(arma::ones<arma::vec>(nb_level));
  //arma::vec up_gauss(nb_level);
  //arma::vec dw_gauss(nb_level);
  if (compute_v == "full" || compute_v == "diag"){
    if (compute_v == "full"){
      //V = compute_full_V(signal_modwt) // missing in action. Roberto's code I think had it.
    } 
    if (compute_v == "diag"){
      
      unsigned int num_field = signal_modwt.n_elem;
      
      arma::vec Aj(signal_modwt.n_elem);
      
      for(unsigned int i = 0; i < num_field; i++){
        // Autocovariance using the Discrete Fourier Transform
        arma::vec temp = dft_acf(signal_modwt(i));
        
        // Sum(V*V) - first_element^2 /2
        Aj(i) = dot(temp,temp) - temp(0)*temp(0)/2;
      }
      // Create diagnoal matrix (2 * Aj / length(modwt_d1)). Note: All modwt lengths are the same. Update if dwt is used.      
      V(0) = diagmat(2 * Aj / signal_modwt(0).n_elem);
      
      //if(robust){
      V(1) = 1/eff*V(0);
      //}
    }
    
    // Compute confidence intervals
    //up_gauss = vmod.col(0) + R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
    //dw_gauss = vmod.col(0) - R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
  }
  
  return V;
}

//' @title Computes the (MODWT) wavelet covariance matrix using Chi-square confidence interval bounds
//' @description Calculates the (MODWT) wavelet covariance matrix using Chi-square confidence interval bounds
//' @param ci_hi A \code{vec} that contains the upper confidence interval points.
//' @param ci_lo A \code{vec} that contains the lower confidence interval points.
//' @return A diagonal matrix.
//' @keywords internal
//' @examples
//' \dontrun{
//' x=runif(100)
//' y=x+3
//' fast_cov_cpp(y,x)
//' }
// [[Rcpp::export]]
arma::mat fast_cov_cpp(const arma::vec& ci_hi, const arma::vec& ci_lo){
  return arma::diagmat(arma::square(ci_hi-ci_lo));
}