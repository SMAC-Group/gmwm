#include <RcppArmadillo.h>
#include "wave_variance.h"
// Uses filters
#include "wv_filters.h"

// Uses brick_wall
#include "dwt.h"

// Uses robust components...
#include "robust_components.h"

using namespace Rcpp;


//' @title Generate eta3 confidence interval
//' @description Computes the eta3 CI
//' @param y A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param dims A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @examples
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, "periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", p = 0.025)
// [[Rcpp::export]]
arma::mat ci_eta3(arma::vec y,  arma::vec dims, double p) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);

    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2,i+1),1.0);
      out(i,1) = eta3 * y(i)/R::qchisq(1-p, eta3, 1, 0); // Lower CI
      out(i,2) = eta3 * y(i)/R::qchisq(p, eta3, 1, 0); // Upper CI
    }

    out.col(0) = y;

    return out;
}

//' @title Generate eta3 robust confidence interval
//' @description Computes the eta3 robust CI
//' @param y A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param dims A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Robust Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @examples
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, boundary="periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw, robust = TRUE,  eff = 0.6)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", p = 0.025, robust = TRUE, eff = 0.6)
// [[Rcpp::export]]
arma::mat ci_eta3_robust(arma::vec y, arma::vec dims, double p, double eff) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);
    eff = sqrt(eff);
    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2,i+1),1.0);
      out(i,1) = eff * eta3 * y(i)/(R::qchisq(1-p, eta3, 1, 0)); // Lower CI
      out(i,2) = eta3 * y(i)/(eff*R::qchisq(p, eta3, 1, 0)); // Upper CI
    }

    out.col(0) = y;

    return out;
}

//' @title Generate a Confidence intervval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a chi-squared confidence interval
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param y A \code{vec} that contains the wave variance.
//' @param type A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level.
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @details 
//' This function can be expanded to allow for other confidence interval calculations.
//' @examples
//' set.seed(1337)
//' x = rnorm(100)
//' # Uses the internal MODWT function not associated with an S3 class.
//' decomp = modwt_cpp(x, "haar", 4, boundary="periodic")
//' signal_modwt_bw = brick_wall(decomp, select_filter("haar"), "modwt")
//' y = wave_variance(signal_modwt_bw)
//' ci_wave_variance(signal_modwt_bw, y, type = "eta3", p = 0.025)
// [[Rcpp::export]]
arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& y,
                            std::string type = "eta3", double p = 0.025, bool robust = false, double eff = 0.6){
    
  unsigned int nb_level = y.n_elem;
  arma::vec dims(nb_level);
  
  for(unsigned int i = 0; i < nb_level; i++){
    dims(i) = signal_modwt_bw(i).n_elem;
  }

  arma::mat out(y.n_elem , 3);

  
  if(type == "eta3"){
      if(robust){
        out = ci_eta3_robust(y, dims, p, eff);
      }
      else{
        out = ci_eta3(y, dims,p);  
      }
  }
  else{
      stop("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}


//' @title Generate a Wave Variance for a Univariate Time Series
//' @description Computes an estimate of the wave variance
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{vec} that contains the wave variance.
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
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @param ci_type A \code{String} indicating the confidence interval being calculated. Valid value: "eta3"
//' @param strWavelet A \code{String} indicating the type of wave filter to be applied. Must be "haar"
//' @return A \code{mat} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @details 
//' This function powers the wvar object. It is also extendable...
//' @examples
//' x=rnorm(100)
//' decomp = modwt(x)
//' wvar_cpp(decomp$data, robust = FALSE)
// [[Rcpp::export]]
arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt,
                   bool robust=false, double eff=0.6, double p = 0.025, 
                   std::string ci_type="eta3", std::string strWavelet="haar") {
  
  // MODWT transform
  arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, select_filter(strWavelet), "modwt");
  
  // Wavelet Variance
  arma::vec y = wave_variance(signal_modwt_bw, robust, eff);
  
  // Confidence Interval
  return ci_wave_variance(signal_modwt_bw, y, ci_type, p, robust, eff);
}

//' @title Computes the MODWT scales
//' @description Calculates the MODWT scales
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @return A \code{vec} that contains 2^1, ... , 2^J
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