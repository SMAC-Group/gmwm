#include <RcppArmadillo.h>

//#include "summary_inf_mod.h"

// Inference goodies
#include "inference.h"

// Model Selection criterion
#include "analytical_matrix_derivatives.h"

// Model Selection criterion
#include "model_selection.h"

using namespace Rcpp;


// [[Rcpp::export]]

arma::field<arma::mat> get_summary(const arma::vec& theta,
                                   const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                   std::string model_type, 
                                   const arma::vec& wv_empir, const arma::vec& theo, const arma::vec& scales,
                                   const arma::mat& V, const arma::mat& omega,
                                   unsigned int N, double alpha, 
                                   bool inference, bool modelselect, bool ci_bootstrap, unsigned int B, bool robust, double eff){
  // Inference test result storage
  arma::vec gof_test;
  arma::mat ci_inf; 
  arma::vec score;
  
  // Take derivatives
  arma::mat D = derivative_first_matrix(theta, desc, objdesc, scales);
  
  // Generate inference information
  if(inference){
    // Obtain a confidence interval for the parameter estimates AND calculate chisq goodness of fit
    arma::field<arma::mat> cat = inference_summary(theta, 
                                                   desc,
                                                   objdesc,
                                                   model_type,
                                                   scales,
                                                   D, 
                                                   V, omega,
                                                   wv_empir, N, alpha , ci_bootstrap, B, robust, eff);
    
    ci_inf = cat(0);
    gof_test = cat(1);
  }
  
  if(modelselect){
    
    // Take derivatives
    arma::mat At_j = derivative_second_matrix(theta, desc, objdesc, scales);
    
    // Obtain the difference
    arma::vec diff = theo - wv_empir;
    
    // Calculate the model score according to model selection criteria paper
    score = model_score(D, At_j, omega, V,  diff, N);
  }

  // Export information back
  arma::field<arma::mat> out(3);
  out(0) = ci_inf;
  out(1) = gof_test;
  out(2) = score;
  
  return out;
  
}