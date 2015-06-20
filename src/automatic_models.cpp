#include <RcppArmadillo.h>

#include "bootstrappers.h"
#include "gmwm_logic.h"

// Used for scales_cpp
#include "wave_variance.h"

#include "analytical_matrix_derivatives.h"

#include "model_selection.h"

#include "ts_model_cpp.h"

//#include "automatic_models.h"
using namespace Rcpp;

//[[Rcpp::export]]
std::map<int, std::vector<std::string> > master_model(){
  
  std::map<int, std::vector<std::string> > models;
  
  // 4AR1() + WN()
  models[0].push_back("AR1");
  models[0].push_back("AR1");
  models[0].push_back("AR1");
  models[0].push_back("AR1");
  models[0].push_back("WN");
  
  // 4AR1()
  models[1].push_back("AR1");
  models[1].push_back("AR1");
  models[1].push_back("AR1");
  models[1].push_back("AR1");
  
   models[2].push_back("AR1");
  models[2].push_back("AR1");
  models[2].push_back("AR1");
  models[2].push_back("WN");
  
  // 3AR1()
  models[3].push_back("AR1");
  models[3].push_back("AR1");
  models[3].push_back("AR1");
  
  
  // 2AR1()+WN()
  models[4].push_back("AR1");
  models[4].push_back("AR1");
  models[4].push_back("WN");
  
  // 2AR1()
  models[5].push_back("AR1");
  models[5].push_back("AR1");
  
  // AR1()
  models[6].push_back("AR1");
  
  // WN()
  models[7].push_back("WN");

  return models;
}

// [[Rcpp::export]]
arma::mat auto_select(const arma::vec& data, 
                      std::string model_type,
                      double alpha, 
                      std::string compute_v, unsigned int K, unsigned int H,
                      unsigned int G, 
                      bool robust, double eff, bool bs_optimism){
  
  // Number of data points
  unsigned int N = data.n_elem;
  
  // Keep track of models parsed
  unsigned int count = 0;
  
  // Create a map
  std::map<int, std::vector<std::string> > models = master_model();
  
  // Make an iterator to iterator through it
  std::map<int, std::vector<std::string> >::iterator iter;
  iter = models.begin(); 
  
  // Get the first model
  std::vector<std::string> desc = iter->second;
  
  // Increment iterator
  iter++;
  
  // Build the fields off of the first model's description
  arma::vec theta = model_theta(desc);
  arma::field<arma::vec> objdesc = model_objdesc(desc); 
  
  // Build matrix to store results
  arma::mat results(2, models.size());
  
  // Obtain the largest models information
  arma::field<arma::mat> master = gmwm_master_cpp(data, 
                                                theta,
                                                desc,
                                                objdesc, 
                                                model_type, 
                                                true, //starting
                                                alpha, 
                                                "fast", // compute V
                                                K, H,
                                                G, 
                                                robust, eff);
    

  // Theta update
  theta = master(0);
  
  // Define WV Empirical
  arma::vec wv_empir = master(2);
  
  // Get bootstrapped V
  arma::mat V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false); // Bootstrapped V (largest model)
  
  // Get the original "FAST" matrix
  arma::mat orgV = master(6); // Original V
  
  // Take the inverse
  arma::mat omega = inv(orgV); // Original V => Omega
  
  // Get expect_diff for guesses with DR
  double expect_diff = arma::as_scalar(master(7));
  
  // Obtain the theoretical WV
  arma::vec theo = master(8);
  
  // Obtain the obj_value of the function
  double obj_value = arma::as_scalar(master(10));
  

  
  // Calculate the values of the Scales 
  arma::vec scales = scales_cpp(floor(log2(N)));
  

  
  /* A note to someone in the future...
   * Yes, there is a difference in order between the diff (wv_empir-theo) for D_matrix
   *  and the model_score diff (theo-wv_empir).
   */
  if(bs_optimism){
    
    arma::mat cov_nu_nu_theta = optimism_bootstrapper(theta,
                                                      desc, objdesc,
                                                      scales, model_type, 
                                                      N, robust, eff, alpha,
                                                      H);
    
    double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
    
    arma::vec temp(2);
    temp(0) = obj_value + optimism;
    temp(1) = optimism;
    results.col(count) = temp;
  }else{
    // Take derivatives
    arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);
    
    // Create the D Matrix (note this is in the analytical_matrix_derivaties.cpp file)
    arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(wv_empir - theo));
    
    // Calculate the model score according to model selection criteria paper
    results.col(count) = model_score(A, D, omega, V,  obj_value);
  }
  
  count++;
  
  while(iter != models.end()){
    // Get the first model
    std::vector<std::string> desc = iter->second;
    
    
    // Build the fields off of the first model's description
    theta = model_theta(desc);
    objdesc = model_objdesc(desc); 
    
    // Run the update version of the GMWM
    arma::field<arma::mat> update = gmwm_update_cpp(theta,
                                                      desc, 
                                                      objdesc, 
                                                      model_type, 
                                                      N, expect_diff, 
                                                      orgV, scales, wv_empir,
                                                      true, //starting
                                                      "fast", 
                                                      K,
                                                      H,
                                                      G, 
                                                      robust,
                                                      eff);
    
    
    // Theta update
    theta = update(0);
    
    // Update theo
    theo = update(3);
    
    // Update objective function
    obj_value = arma::as_scalar(update(5));
    
    /* A note to someone in the future...
     * Yes, there is a difference in order between the diff (wv_empir-theo) for D_matrix
     *  and the model_score diff (theo-wv_empir).
     */    
    if(bs_optimism){
      
      arma::mat cov_nu_nu_theta = optimism_bootstrapper(theta,
                                                        desc, objdesc,
                                                        scales, model_type, 
                                                        N, robust, eff, alpha,
                                                        H);
      
      double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
      
      arma::vec temp(2);
      temp(0) = obj_value + optimism;
      temp(1) = optimism;
      results.col(count) = temp;
    }else{
      // Take derivatives
      arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);
      
      // Create the D Matrix (note this is in the analytical_matrix_derivaties.cpp file)
      arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(wv_empir - theo));
      
      // Calculate the model score according to model selection criteria paper
      results.col(count) = model_score(A, D, omega, V,  obj_value);
    }
    
    // Increment iterator
    iter++;
    
    // Increase count
    count++;
  }
  
  return results;
}