#include <RcppArmadillo.h>

//#include "summary_inf_mod.h"

// Inference goodies
#include "inference.h"

// Model Selection criterion
#include "analytical_matrix_derivatives.h"

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
                                   bool inference, bool fullV,
                                   bool bs_gof,  bool bs_gof_p_ci, bool bs_theta_est, bool bs_ci, 
                                   unsigned int B){
  
  // Inference test result storage
  arma::vec gof;
  
  arma::mat ci;
  

  // Bootstrap results
  arma::vec sd;

  arma::vec bs_obj_values;
    

  // Determine the type of bootstrapper needed.
  
  if(!(bs_ci || bs_gof) & !fullV & inference){
     V = cov_bootstrapper(theta,
                       desc, objdesc,
                       N, robust, eff,
                       B, true);
  }else if(bs_ci || bs_gof){
    arma::field<arma::mat> bs = all_bootstrapper(theta,
                                                 desc, objdesc,
                                                 scales, model_type, 
                                                 N, robust, eff, alpha, B);
    
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
    }else{
      arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);
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


  // Export information back
  arma::field<arma::mat> out(2);
  out(0) = ci;
  out(1) = gof;
  
  return out;
  
}