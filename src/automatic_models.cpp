#include <RcppArmadillo.h>

#include "bootstrappers.h"
#include "inference.h"
#include "gmwm_logic.h"

// Used for scales_cpp
#include "wave_variance.h"

#include "analytical_matrix_derivatives.h"

#include "model_selection.h"

#include "ts_model_cpp.h"

//#include "automatic_models.h"
using namespace Rcpp;


// ---- START helper functions

//' @title Build List of Unique Models
//' @description Creates a set containing unique strings. 
//' @param combs A \code{mat} that is a binary matrix (0,1) containing the combinations of different variables.
//' @param x A \code{vec<string>} that contains a list of model descriptors.
//' @return A \code{set<string>} that contains the list of unique models.
// [[Rcpp::export]]
std::set<std::vector<std::string > > build_model_set(const arma::mat& combs, std::vector <std::string> x) {
  
  std::set<std::vector<std::string > > models;
  for(unsigned int i = 0; i < combs.n_rows; i++) {
    std::vector< std::string > tmp;
    for(unsigned int j = 0; j < combs.n_cols; j++){
      if(combs(i,j) ==  1){ 
        tmp.push_back( x[j] );
      }
    }
    models.insert(tmp);
  }

  return models;
}


//' @title Find the Common Denominator of the Models
//' @description Determines the common denominator among models
//' @param x A \code{vector< vector<string> >} that contains all possible models under consideration
//' @return A \code{vector<string>} that contains the terms of the common denominator of all models
// [[Rcpp::export]]
std::vector<std::string> find_full_model(std::vector<std::vector<std::string> > x ){
  
  // Set up iterators to go over the sets
  std::vector<std::vector<std::string> >::const_iterator it;
  std::vector<std::string>::const_iterator it2;
  
  
  // Figure out all common terms
  
  // AR1s can have an infinite amount of combinations
  unsigned int maxAR1s = 0;
  
  // In the mean time, WN, RW, QN, and DR, can only appear once in a model. 
  bool WN = false;
  bool RW = false;
  bool QN = false;
  bool DR = false;
  
  
  // Begin iterating through the set. 
  for(it = x.begin(); it != x.end(); ++it)
  {
    
    // Create an internal counter of AR1s for a given vector
    unsigned int num_AR1s = 0;
    
    // Iterate through the vector 
    for (it2 = (*it).begin(); it2 != (*it).end(); ++it2){

      if(*it2 == "AR1"){
        num_AR1s++; // For each AR1 detected, increment by 1. 
        
      }else if(*it2 == "WN"){
        
        if(!WN){
          WN = true; // If WN has yet to be set, set it here. 
        }
        
      }else if(*it2 == "RW"){
        
        if(!RW){
          RW = true;
        }
        
      }else if(*it2 == "QN"){
        
        if(!QN){
          QN = true;
        }
        
      }else if(*it2 == "DR"){
        
        if(!DR){
          DR = true;
        }
        
      }
      // end if block
      
      // Does this model have more AR1s than previous models? Y/N?
      if(num_AR1s > maxAR1s){
        maxAR1s = num_AR1s;
      }
      
    }
    // end inner for 
    
  }
  // end outer for
  
  // Create a vector holding all the model terms
  std::vector<std::string> out(maxAR1s+WN+RW+QN+DR);
  
  
  // Begin assigning the terms
  unsigned int i;
  
  for(i = 0; i < maxAR1s; i++){
    out[i] = "AR1";
  }
  
  if(WN){
    out[i] = "WN";
    i++;
  }
  
  if(RW){
    out[i] = "RW";
    i++;
  }
  
  if(QN){
    out[i] = "QN";
    i++;
  }
  
  if(DR){
    out[i] = "DR";
    i++;
  }
  
  return out;
}




arma::rowvec bs_optim_calc(const arma::vec& theta,
                        const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                        std::string model_type, const arma::vec& scales, const arma::mat& omega, unsigned int N,
                        double obj_value, double alpha,
                        std::string compute_v, 
                        unsigned int K, unsigned int H, unsigned int G, 
                        bool robust, double eff){
  
  std::cout << "Calling Opt Bootstrap" << std::endl;
  
  arma::field<arma::mat> bso = opt_n_gof_bootstrapper(theta,
                                                      desc, objdesc,
                                                      scales, model_type, 
                                                      N, robust, eff, alpha,
                                                      H);
  arma::mat cov_nu_nu_theta = bso(0);
  
  arma::mat bs_obj_values = bso(1);
  
  std::cout << "Assigned bootstrapped values" << std::endl;
  
  double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
  
  std::cout << "Creating output vector" << std::endl;
  
  arma::rowvec temp(4);
  
  temp(0) = obj_value;
  
  Rcpp::Rcout << "First spot in temp is: " << temp(0) << std::endl;
  temp(1) = optimism;
  
  temp(2) = obj_value + optimism;
  Rcpp::Rcout << "Second spot in temp is: " << temp(1) << std::endl;
  
  std::cout << "Assigned obj and criterion... trying bootstrap" << std::endl;
  
  temp(3) = arma::as_scalar(bootstrap_gof_test(obj_value, bs_obj_values, alpha, false).row(0));
  
  std::cout << "Bootstrap assigned" << std::endl;
  
  Rcpp::Rcout << "returned temp values" << temp << std::endl;
  
  return temp;
}


arma::rowvec asympt_calc(const arma::vec& theta,
                         const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                         std::string model_type, const arma::vec& scales, const arma::mat& V, const arma::mat& omega,
                         const arma::vec& wv_empir, const arma::vec& theo, double obj_value){
  
  // Take derivatives
  arma::mat A = derivative_first_matrix(theta, desc, objdesc, scales);

  /* A note to someone in the future...
   * Yes, there is a difference in order between the diff (wv_empir-theo) for D_matrix
   *  and the model_score diff (theo-wv_empir).
   */
  
  // Create the D Matrix (note this is in the analytical_matrix_derivaties.cpp file)
  arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(wv_empir - theo));
  
  arma::rowvec temp(4);
  
  arma::vec result = model_score(A, D, omega, V,  obj_value);
  
  temp(0) = obj_value;
  temp(1) = arma::as_scalar(result.row(1));
  temp(2) = arma::as_scalar(result.row(0));
  temp(3) = arma::as_scalar( gof_test(theta, desc, objdesc, model_type, scales, V, wv_empir).row(1) ); // GoF
  
  return temp;
}


// ---- End helper functions

arma::field<arma::field<arma::mat> > model_select(const arma::mat& data,
                          const std::set<std::vector<std::string > >& models,
                          const std::vector< std::string >& full_model,
                          std::string model_type,
                          bool bs_optimism,
                          double alpha,
                          std::string compute_v, 
                          unsigned int K, unsigned int H, unsigned int G, 
                          bool robust, double eff){
  
  // Number of data points
  unsigned int N = data.n_rows;
  
  // Number of models
  unsigned int num_models = models.size(); 
  
  // Make an iterator to iterator through it
  std::set<std::vector<std::string > > ::const_iterator iter;
  iter = models.begin(); 
  
  // Get the first model
  std::vector<std::string> desc = full_model;
  
  std::cout << "We are currently running the full model" << std::endl;
  for(unsigned int i = 0; i<desc.size(); i++){ std::cout << desc[i] << " "; }
  std::cout << std::endl << "End full model" << std::endl;
    
  // Find where the results should be input. (No protection needed, we know it is in the matrix)
  unsigned int full_model_index = std::distance(models.begin(),models.find(full_model));
  std::cout << "The full model location is " << full_model_index << std::endl;
  
  
  // Build the fields off of the first model's description
  arma::vec theta = model_theta(desc);
  arma::field<arma::vec> objdesc = model_objdesc(desc); 
  
  // Build matrix to store results
  arma::mat results(num_models, 4);
  
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
  
  Rcpp::Rcout << "The value of the GMWM object is " << master << std::endl;
  
  // Theta update
  theta = master(0);
  
  // Define WV Empirical
  arma::vec wv_empir = master(2);
  
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
  
  std::cout << "Objective function for full model: " << obj_value << std::endl;
  
  // Calculate the values of the Scales 
  arma::vec scales = scales_cpp(floor(log2(N)));
  
  // ------------------------------------
  
  // Here we set up specifics that are used not in a specific mode.
  
  // Asymptotic ----
  
  // Get bootstrapped V
  arma::mat V;
  
  // Bootstrap ----
  // Hold optimism result
  arma::mat cov_nu_nu_theta;
  
  // Hold the bootstrapped obj values
  arma::vec bs_obj_values;
  
  if(bs_optimism){
    std::cout << "Testing Optimism bootstrapper" << std::endl;
    results.row(full_model_index) = bs_optim_calc(theta,  desc,  objdesc, model_type, scales, omega, N,
                obj_value, alpha, compute_v, K, H, G, robust, eff);
    std::cout << "Passed Optimism bootstrapper" << std::endl;
  }else{
    
    std::cout << "Calculating the CoV, V Matrix, bootstrapper" << std::endl;
    
    V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false) * sqrt(N); // Bootstrapped V (largest model)
    
    std::cout << "End calculation for the CoV, V Matrix, bootstrapper" << std::endl;
    
    Rcpp::Rcout << "V bootstrapped matrix" << std::endl << V << std::endl;
    
    std::cout << "Calculating the asymptotic model" << std::endl;
    
    // Calculate the model score according to model selection criteria paper
    results.row(full_model_index) = asympt_calc(theta, desc, objdesc, model_type, scales, V, omega, wv_empir, theo, obj_value);
    
    std::cout << "End calculation for the asymptotic model" << std::endl;
  }
  
  // Initialize counter to keep track of values
  unsigned int count = 0;
  
  while(iter != models.end()){
    
    if(full_model_index != count){
      // Get the first model
      desc = *iter;
      
      std::cout << "We are currently running a nested model" << std::endl;
      for(unsigned int i = 0; i<desc.size(); i++){ std::cout << desc[i] << " "; }
      std::cout << std::endl << "End nested model" << std::endl;
      
      // Build the fields off of the first model's description
      theta = model_theta(desc);
      objdesc = model_objdesc(desc); 
      
      // Run the update version of the GMWM
      arma::field<arma::mat> update = gmwm_update_cpp(theta, desc, objdesc, model_type, 
                                                      N, expect_diff, 
                                                      orgV, scales, wv_empir,
                                                      true, //starting
                                                      "fast", 
                                                      K,H,G, 
                                                      robust,eff);
      
      
      // Theta update
      theta = update(0);
      
      // Update theo
      theo = update(3);
      
      // Update objective function
      obj_value = arma::as_scalar(update(5));
      
      std::cout << "Objective function for nested model: " << obj_value << std::endl;
      
      if(bs_optimism){
        results.row(count) = bs_optim_calc(theta,  desc,  objdesc, model_type, scales, omega, N,
                    obj_value, alpha, compute_v, K, H, G, robust, eff);
      }else{
        // Calculate the model score according to model selection criteria paper
        results.row(count) = asympt_calc(theta, desc, objdesc, model_type, scales, V, omega, wv_empir, theo, obj_value);
      }
    }
    // end if
    
    // Increment iterator
    iter++;
    
    // Increase count
    count++;
  }
  
  std::cout << "Finished computing results!" << std::endl;
  Rcpp::Rcout << results << std::endl;
  arma::field< arma::field<arma::mat> > out(1);
  
  arma::field<arma::mat> ms(1);
  ms(0) = results;
  
  out(0) = ms;
  return out;
}



// [[Rcpp::export]]
arma::field< arma::field<arma::field<arma::mat> > >  auto_imu(const arma::mat& data,
                                              const arma::mat& combs,
                                              const std::vector< std::string >&  full_model,
                                              double alpha, 
                                              std::string compute_v, std::string model_type, 
                                              unsigned int K, unsigned int H, unsigned int G, 
                                              bool robust, double eff, bool bs_optimism){
  
  
  
  
  // Number of columns to process
  unsigned int V = data.n_cols;
  
  // Create a set of unique models.
  std::set<std::vector<std::string > > models = build_model_set(combs, full_model);
  
  arma::field< arma::field<arma::field<arma::mat> > > h(V);
  
  for(unsigned int i = 0; i < V; i++){
    h(i) = model_select(data.col(i),
      models,
      full_model,
      model_type,
      bs_optimism,
      alpha,
      compute_v, 
      K, H, G, 
      robust, eff);
  }
  
  return h;
}
