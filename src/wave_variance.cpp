/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Q Public License included within the packages source
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Q Public License
 * along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#include <RcppArmadillo.h>
#include "wave_variance.h"
// Uses filters
#include "wv_filters.h"

// Uses brick_wall
#include "dwt.h"

// Uses robust components...
#include "robust_components.h"

//' @title Generate eta3 confidence interval
//' @description Computes the eta3 CI
//' @param y A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param dims A \code{String} indicating the confidence interval being calculated.
//' @param alpha_ov_2 A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @keywords internal
//' @examples
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, "periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", alpha_ov_2 = 0.025)
// [[Rcpp::export]]
arma::mat ci_eta3(arma::vec y,  arma::vec dims, double alpha_ov_2) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);

    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2,i+1),1.0);
      out(i,1) = eta3 * y(i)/R::qchisq(1-alpha_ov_2, eta3, 1, 0); // Lower CI
      out(i,2) = eta3 * y(i)/R::qchisq(alpha_ov_2, eta3, 1, 0); // Upper CI
    }

    out.col(0) = y;

    return out;
}

//' @title Generate eta3 robust confidence interval
//' @description Computes the eta3 robust CI
//' @param wv_robust A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param wv_ci_class A \code{mat} that contains the CI mean, CI Lower, and CI Upper
//' @param alpha_ov_2 A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Robust Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @details
//' Within this function we are scaling the classical 
//' @keywords internal
//' @examples
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, boundary="periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw, robust = TRUE,  eff = 0.6)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", alpha_ov_2 = 0.025, robust = TRUE, eff = 0.6)
// [[Rcpp::export]]
arma::mat ci_eta3_robust(arma::vec wv_robust, arma::mat wv_ci_class, double alpha_ov_2, double eff) {
    unsigned int num_elem = wv_robust.n_elem;

    arma::mat out(num_elem, 3);
    
    double q1 = R::qnorm(1-alpha_ov_2, 0.0, 1.0, true, false);
    
    double coef = ((-1.0*q1-q1) * sqrt(1.0/eff)) / (-1.0*q1-q1);

    for(unsigned int i = 0; i<num_elem;i++){
      
      double wv_ci = wv_ci_class(i,0); // store WV for i-th case
      
      
      double lci = (wv_ci - wv_ci_class(i,1))/wv_ci; // lower classical ci scale
      double uci = (wv_ci_class(i,2) - wv_ci)/wv_ci; // upper classical ci scale
      
      
      double wv_ri = wv_robust(i); // store WV for i-th case
      
      // Avoid negative result
      lci = wv_ri - lci*coef*wv_ri;
      if(lci > 0){
        out(i,1) = lci;        
      }else{
        out(i,1) = wv_ci_class(i,1)/2; // mean of (low_ci, 0) /2
      }

      out(i,2) = wv_ri + uci*coef*wv_ri;
    }

    out.col(0) = wv_robust;

    return out;
}

//' @title Generate a Confidence intervval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a chi-squared confidence interval
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param wv A \code{vec} that contains the wave variance.
//' @param type A \code{String} indicating the confidence interval being calculated.
//' @param alpha_ov_2 A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level.
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @keywords internal
//' @details 
//' This function can be expanded to allow for other confidence interval calculations.
//' @examples
//' set.seed(1337)
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, boundary="periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", alpha_ov_2 = 0.025)
// [[Rcpp::export]]
arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& wv,
                            std::string type = "eta3", double alpha_ov_2 = 0.025, bool robust = false, double eff = 0.6){
    
  unsigned int nb_level = wv.n_elem;
  arma::vec dims(nb_level);
  
  for(unsigned int i = 0; i < nb_level; i++){
    dims(i) = signal_modwt_bw(i).n_elem;
  }

  arma::mat out(nb_level , 3);

  
  if(type == "eta3"){
      
      if(!robust){
        out = ci_eta3(wv, dims, alpha_ov_2);  
      }else{
        // per the WV Robust change... 
        // We need to obtain the classical CI first, then modify it.
        arma::vec wv_class = wave_variance(signal_modwt_bw, false, eff); // Requires the next function....
        arma::mat wv_ci_class = ci_eta3(wv_class, dims, alpha_ov_2);  // calculate the CI
    
        // wv is the wave robust
        out = ci_eta3_robust(wv, wv_ci_class, alpha_ov_2, eff);
      }
  }
  else{
    Rcpp::stop("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}


//' @title Generate a Wave Variance for a Univariate Time Series
//' @description Computes an estimate of the wave variance
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{vec} that contains the wave variance.
//' @keywords internal
//' @examples
//' set.seed(1337)
//' x = rnorm(100)
//' signal_modwt_bw = brick_wall(modwt_cpp(x), haar_filter())
//' wave_variance(signal_modwt_bw)
//' 
//' wave_variance(signal_modwt_bw, robust = TRUE, eff = 0.6)
// [[Rcpp::export]]
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw, bool robust = false, double eff = 0.6){
  
  unsigned int nb_level = signal_modwt_bw.n_elem;
  arma::vec y(nb_level);
  
  if(robust){
    // Robust wavelet variance estimation
    for(unsigned int i=0; i < nb_level; i++){
      arma::vec wav_coef = sort(signal_modwt_bw(i));
      y(i) = sig_rob_bw(wav_coef, eff);
    }
  }else{
    // Classical wavelet variance estimation
    for(unsigned int i=0; i < nb_level;i++){
      arma::vec temp = signal_modwt_bw(i);
      y(i) = dot(temp,temp)/temp.n_elem;
    }
  }
  
  return y;
}


//' @title Computes the (MODWT) wavelet variance
//' @description Calculates the (MODWT) wavelet variance
//' @param signal_modwt A \code{field<vec>} that contains the modwt decomposition.
//' @param robust A \code{boolean} that triggers the use of the robust estimate.
//' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
//' @param alpha A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @param ci_type A \code{String} indicating the confidence interval being calculated. Valid value: "eta3"
//' @param strWavelet A \code{String} indicating the type of wave filter to be applied. Must be "haar"
//' @return A \code{mat} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @keywords internal
//' @details 
//' This function powers the wvar object. It is also extendable...
//' @examples
//' x=rnorm(100)
//' decomp = modwt(x)
//' wvar_cpp(decomp$data, robust = FALSE)
// [[Rcpp::export]]
arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt,
                   bool robust=false, double eff=0.6, double alpha = 0.05, 
                   std::string ci_type="eta3", std::string strWavelet="haar") {
  double alpha_ov_2 = alpha/2.0;
  // MODWT transform
  arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, select_filter(strWavelet), "modwt");
  
  // Wavelet Variance
  arma::vec y = wave_variance(signal_modwt_bw, robust, eff);
  
  // Confidence Interval
  return ci_wave_variance(signal_modwt_bw, y, ci_type, alpha_ov_2, robust, eff);
}

//' @title Computes the MODWT scales
//' @description Calculates the MODWT scales
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @return A \code{vec} that contains 2^1, ... , 2^J
//' @keywords internal
//' @details 
//' Used in wvar object.
//' @examples
//' scales_cpp(5)
// [[Rcpp::export]]
arma::vec scales_cpp(unsigned int nb_level){
  // Define scales
  arma::vec scales(nb_level);
  for(unsigned int i=0; i< nb_level;i++){
    scales(i) = pow(2,i+1);
  }
  return scales;
}