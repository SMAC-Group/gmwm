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


//' @title Conversion function of Vector to Set
//' @description Converts a vector into a set
//' @param x A \code{vec<vec<string>>} that contains a list of model descriptors.
//' @return A \code{set<vector<string>>} that contains the list of unique models.
// [[Rcpp::export]]
std::set<std::vector<std::string> > vector_to_set(std::vector<std::vector<std::string > > model_str){
  
  // Number of columns to process
  // Create a set of unique models.
  
  std::set<std::vector<std::string > > models;
  for(std::vector<std::vector<std::string> >::const_iterator it = model_str.begin(); it != model_str.end(); ++it) {
    models.insert(*it);
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
  
  arma::field<arma::mat> bso = opt_n_gof_bootstrapper(theta,
                                                      desc, objdesc,
                                                      scales, model_type, 
                                                      N, robust, eff, alpha,
                                                      H);
  arma::mat cov_nu_nu_theta = bso(0);
  
  arma::mat bs_obj_values = bso(1);
  

  double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
  

  arma::rowvec temp(4);
  
  temp(0) = obj_value;
  
  temp(1) = optimism;
  
  temp(2) = obj_value + optimism;

  temp(3) = arma::as_scalar(bootstrap_gof_test(obj_value, bs_obj_values, alpha, false).row(0));

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
  arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(theo - wv_empir));
  
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
  
  // Find where the results should be input. (No protection needed, we know it is in the matrix)
  unsigned int full_model_index = std::distance(models.begin(),models.find(full_model));
  
  // Build the fields off of the first model's description
  arma::vec theta = model_theta(desc);
  arma::field<arma::vec> objdesc = model_objdesc(desc); 
  
  // Build matrix to store results
  arma::mat results(num_models, 4);
  
  std::cout << "Processing model 1 out of " << num_models << std::endl;
  
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
    results.row(full_model_index) = bs_optim_calc(theta,  desc,  objdesc, model_type, scales, omega, N,
                obj_value, alpha, compute_v, K, H, G, robust, eff);
  }else{
    
    std::cout << "One time computation of the V matrix via bootstrap... Please stand by." << std::endl;
    V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false); // Bootstrapped V (largest model)
    
    // Calculate the model score according to model selection criteria paper
    results.row(full_model_index) = asympt_calc(theta, desc, objdesc, model_type, scales, V, omega, wv_empir, theo, obj_value);
    
  }
  
  // Initialize counter to keep track of values
  unsigned int count = 0;
  
  unsigned int countModels = 1;
  
  while(iter != models.end()){
    
    if(full_model_index != count){
      countModels++;
      std::cout << "Processing model " << countModels << " out of " << num_models << std::endl;
      // Get the first model
      desc = *iter;
      
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
  
  arma::field< arma::field<arma::mat> > out(1);
  
  arma::field<arma::mat> ms(1);
  ms(0) = results;
  
  out(0) = ms;
  return out;
}



// [[Rcpp::export]]
arma::field< arma::field<arma::field<arma::mat> > >  rank_models(const arma::vec& data,
                                                                const std::vector<std::vector < std::string > >& model_str, 
                                                                const std::vector< std::string >&  full_model,
                                                                double alpha, 
                                                                std::string compute_v, std::string model_type, 
                                                                unsigned int K, unsigned int H, unsigned int G, 
                                                                bool robust, double eff, bool bs_optimism){
  
  
  std::set<std::vector < std::string > > models = vector_to_set(model_str);
  
  arma::field< arma::field<arma::field<arma::mat> > > h(1);
  
  h(0) = model_select(data,
      models,
      full_model,
      model_type,
      bs_optimism,
      alpha,
      compute_v, 
      K, H, G, 
      robust, eff);
  
  return h;
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
