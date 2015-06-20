#include <RcppArmadillo.h>

//#include "summary_inf_mod.h"

// Inference goodies
#include "inference.h"

// Model Selection criterion
#include "analytical_matrix_derivatives.h"

// Model Selection criterion
#include "model_selection.h"

// Bootstraps!
#include "bootstrappers.h"

using namespace Rcpp;


// [[Rcpp::export]]

arma::field<arma::mat> get_summary(arma::vec theta,
                                   const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                   std::string model_type, 
                                   const arma::vec& wv_empir, const arma::vec& theo, const arma::vec& scales,
                                   arma::mat V, const arma::mat& omega, double obj_value,
                                   unsigned int N, double alpha,
                                   bool robust, double eff, 
                                   bool inference, bool model_select, bool fullV,
                                   bool bs_gof,  bool bs_gof_p_ci, bool bs_theta_est, bool bs_ci, bool bs_optimism, 
                                   unsigned int B){
  
  // Inference test result storage
  arma::vec gof;
  
  arma::mat ci;
  
  arma::vec score;
  

  // Bootstrap results
  arma::vec sd;
  
  arma::mat cov_nu_nu_theta;

  arma::vec bs_obj_values;
    
    // Take derivatives
  arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);
  
  
  if(!(bs_ci || bs_optimism || bs_gof) & !fullV & (inference || model_select)){
    Rcout << "stand alone" << std::endl;
     V = cov_bootstrapper(theta,
                       desc, objdesc,
                       N, robust, eff,
                       B, true);
  }else if(bs_ci || bs_optimism || bs_gof){
    Rcout << "BS model" << std::endl;
    
    arma::field<arma::mat> bs = all_bootstrapper(theta,
                                                 desc, objdesc,
                                                 scales, model_type, 
                                                 N, robust, eff, alpha, B);
    
    if(bs_optimism){
      cov_nu_nu_theta = bs(0);
    }

    if(!fullV){
      V = bs(1);
    }
    
    if(bs_theta_est){
      theta = bs(2);
    }
    
    if(bs_ci){
      sd = bs(3);
    }
    
    if(bs_gof){
      bs_obj_values = bs(4);
    }
  }

  // Generate inference information
  if(inference){
    // Obtain a confidence interval for the parameter estimates AND calculate chisq goodness of fit
    if(bs_ci){
      ci = format_ci(theta, sd, alpha);
    }else if(!bs_ci){
      ci = theta_ci(theta, A, V, omega, alpha);
    }
    
    if(bs_gof){
      arma::vec temp(3);
      
      temp(0) = sum(obj_value < bs_obj_values)/double(B);
      if(bs_gof_p_ci){
        temp.rows(1,2) = boot_pval_gof(obj_value, bs_obj_values, 1000, alpha);
      }
      gof = temp; 
    }else{
      gof = gof_test(theta, desc, objdesc, model_type, scales, V, wv_empir);
    }

  }

  if(model_select){
    
    /* A note to someone in the future...
     * Yes, there is a difference in order between the diff (wv_empir-theo) for D_matrix
     *  and the model_score diff (theo-wv_empir).
     */
    if(bs_optimism){
      
      Rcout << "BS" << std::endl;
      arma::vec temp(2);
      
      double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
      
      temp(0) = obj_value + optimism;
      temp(1) = optimism;
      
      score = temp;
    }else{
      Rcout << "Asymptotic" << std::endl;
      
      // Create the D Matrix (note this is in the analytical_matrix_derivaties.cpp file)
      arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(wv_empir - theo));
    
      // Calculate the model score according to model selection criteria paper
      score = model_score(A, D, omega, V,  obj_value);
    }
  }
  

  // Export information back
  arma::field<arma::mat> out(3);
  out(0) = ci;
  out(1) = gof;
  out(2) = score;
  
  return out;
  
}