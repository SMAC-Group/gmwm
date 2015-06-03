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
