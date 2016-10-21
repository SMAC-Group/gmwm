/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
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

#include "bootstrappers.h"
#include "inference.h"
#include "gmwm_logic.h"

// Used for armadillo manipulations
#include "armadillo_manipulations.h"

// Used for scales_cpp
#include "wave_variance.h"

#include "analytical_matrix_derivatives.h"

#include "model_selection.h"

#include "ts_model_cpp.h"

#include "ts_checks.h"

//#include "automatic_models.h"

// ---- START helper functions

unsigned int count_params(const std::vector<std::string>& desc){    
  std::map<std::string, int> w = count_models(desc);
  
  unsigned int params = 0; 
  for (std::map<std::string, int>::iterator it = w.begin(); it!=w.end(); ++it) {		
    if(it->first == "AR1" || it->first == "GM" || it->first == "MA1"){
      params += 2*it->second;
    }else if(it->first != "ARMA11"){
      params += 1;
    }else{
      params += 3;
    }
  }		
  
  return params;		
} 

//' @title Build List of Unique Models
//' @description Creates a set containing unique strings. 
//' @param combs A \code{mat} that is a binary matrix (0,1) containing the combinations of different variables.
//' @param x A \code{vec<string>} that contains a list of model descriptors.
//' @return A \code{set<string>} that contains the list of unique models.
//' @keywords internal
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



//' Set the RNG Seed from within Rcpp
//' 
//' Within Rcpp, one can set the R session seed without triggering
//' the CRAN rng modifier check. 
//' @param seed A \code{unsigned int} that is the seed one wishes to use. 
//' @return A set RNG scope.
//' @keywords internal
//' @examples
//' set.seed(10)
//' x = rnorm(5,0,1)
//' set_seed(10)
//' y = rnorm(5,0,1)
//' all.equal(x,y, check.attributes = FALSE)
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}


//' @title Conversion function of Vector to Set
//' @description Converts a vector into a set
//' @param x A \code{vec<vec<string>>} that contains a list of model descriptors.
//' @return A \code{set<vector<string>>} that contains the list of unique models.
//' @keywords internal
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
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::string> find_full_model(std::vector<std::vector<std::string> > x ){
  
  // Set up iterators to go over the sets
  std::vector<std::vector<std::string> >::const_iterator it;
  std::vector<std::string>::const_iterator it2;
  
  
  // Figure out all common terms
  
  // AR1s can have an infinite amount of combinations
  unsigned int maxAR1s = 0;
  
  // How about MA1s?
  unsigned int maxMA1s = 0;
  
  // And we add in the ability for more than 1 ARMA11().
  unsigned int maxARMA11s = 0;
  
  // String type to describe internal representation (e.g. AR1 vs. GM).
  std::string type = "AR1";
  
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
    unsigned int num_MA1s = 0;
    unsigned int num_ARMA11s = 0;
    
    // Iterate through the vector 
    for (it2 = (*it).begin(); it2 != (*it).end(); ++it2){
      
      if(*it2 == "AR1" || *it2 == "GM"){
        num_AR1s++; // For each AR1 detected, increment by 1. 
        if(*it2 == "GM"){type="GM";}else{type="AR1";}
      }else if(*it2 == "MA1"){
        num_MA1s++; // For each MA1 detected, increment by 1. 
      }else if(*it2 == "ARMA11"){
        num_ARMA11s++; // For each ARMA11 detected, increment by 1. 
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
      
      if(num_MA1s > maxMA1s){
        maxMA1s = num_MA1s;
      }
      
      if(num_ARMA11s > maxARMA11s){
        maxARMA11s = num_ARMA11s;
      }
      
    }
    // end inner for 
    
  }
  // end outer for
  
  // Create a vector holding all the model terms
  std::vector<std::string> out(maxAR1s+maxMA1s+maxARMA11s +WN+RW+QN+DR);
  
  
  // Begin assigning the terms
  unsigned int i;
  
  // Fill vector with ARs first.
  for(i = 0; i < maxAR1s; i++){
    out[i] = type;
  }
  
  for(unsigned int c = 0; c < maxMA1s; c++){
    out[c+i] = "MA1";
  }
  i += maxMA1s;
  
  for(unsigned int c = 0; c < maxARMA11s; c++){
    out[c+i] = "ARMA11";
  }
  i += maxARMA11s;
  
  if(WN){
    out[i] = "WN";
    i++;
  }
  
  if(QN){
    out[i] = "QN";
    i++;
  }
  
  if(RW){
    out[i] = "RW";
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

arma::field<arma::field<arma::mat> > model_select(arma::vec& data,
                                                  const std::set<std::vector<std::string > >& models,
                                                  const std::vector< std::string >& full_model,
                                                  std::string model_type,
                                                  bool bs_optimism,
                                                  double alpha,
                                                  std::string compute_v, 
                                                  unsigned int K, unsigned int H, unsigned int G, 
                                                  bool robust, double eff, unsigned int seed){
  
  // Number of data points
  unsigned int N = data.n_rows;
  
  // Number of models
  unsigned int num_models = models.size(); 
  
  // Store output from models
  arma::field<arma::field<arma::mat> > model_results(num_models);
  
  // Make an iterator to iterator through it
  std::set<std::vector<std::string > > ::const_iterator iter;
  iter = models.begin(); 
  
  // Get the first model
  std::vector<std::string> desc = full_model;
  
  // Find where the results should be input. (No protection needed, we know it is in the matrix)
  unsigned int full_model_index = std::distance(models.begin(), models.find(full_model));
  
  // Build the fields off of the first model's description
  arma::vec theta = model_theta(desc);
  arma::field<arma::vec> objdesc = model_objdesc(desc); 
  
  // Build matrix to store results
  arma::mat results(num_models, 4);
  
  Rcpp::Rcout << "Processing model 1 out of " << num_models << std::endl;
  
  set_seed(seed);
  
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
  
  // Create WV Matrix
  
  arma::mat wv(wv_empir.n_elem,3);
  wv.col(0) = wv_empir;
  wv.col(1) = master(3);
  wv.col(2) = master(4);
  
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
  
  double dr_slope = arma::as_scalar(master(12));
  
  // ------------------------------------
  
  // Store output from default GMWM object
  arma::field<arma::mat> mod_output(13);
  mod_output(0) = theta;
  mod_output(1) = master(1);
  mod_output(2) = wv_empir;
  mod_output(3) = master(3);
  mod_output(4) = master(4);
  mod_output(5) = master(5);
  mod_output(6) = orgV;
  mod_output(7) = expect_diff;
  mod_output(8) = theo;
  mod_output(9) = master(9);
  mod_output(10) = obj_value;
  mod_output(11) = omega;
  mod_output(12) = dr_slope;
  
  
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
    
    Rcpp::Rcout << "Bootstrapping the covariance matrix... Please stand by." << std::endl;
    V = cov_bootstrapper(theta, desc, objdesc, N, robust, eff, H, false); // Bootstrapped V (largest model)
    
    // Calculate the model score according to model selection criteria paper
    results.row(full_model_index) = asympt_calc(theta, desc, objdesc, model_type, scales, V, omega, wv_empir, theo, obj_value);
    
  }
  
  // Custom GMWM Update Obj
  
  arma::field<arma::mat> model_update_info(6);
  model_update_info(0) = master(0); // theta
  model_update_info(1) = master(1); // guessed_theta
  model_update_info(2) = master(5); // V
  model_update_info(3) = master(8); // theo
  model_update_info(4) = master(9); // decomp_theo
  model_update_info(5) = master(10); // objective value
  
  model_results(full_model_index) = model_update_info;
  
  // Initialize counter to keep track of values
  unsigned int count = 0;
  
  unsigned int countModels = 1;
  
  while(iter != models.end()){
    
    if(full_model_index != count){
      countModels++;
      
      // Set guessing seed
      set_seed(seed);
      
      Rcpp::Rcout << "Processing model " << countModels << " out of " << num_models << std::endl;
      // Get the first model
      desc = *iter;
      
      // Build the fields off of the first model's description
      theta = model_theta(desc);
      objdesc = model_objdesc(desc); 
      
      // Run the update version of the GMWM
      arma::field<arma::mat> update = gmwm_update_cpp(theta, desc, objdesc, model_type, 
                                                      N, expect_diff, dr_slope,
                                                      orgV, scales, wv,
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
      
      model_results(count) = update;
      
    }
    // end if
    
    // Increment iterator
    iter++;
    
    // Increase count
    count++;
  }
  
  
  
  // Only run if in asymptotic mode
  if(!bs_optimism){
    
    // std::set<std::vector<std::string > > ::const_iterator iter2;
    // 
    // iter = models.begin(); 
    // 
    // /*
    // * Iterate through all models
    // * If i != j, then 
    // * IF (Complex_i > Complex_j && Obj_i > Obj_j){
    // *  IF(Crit_i < Crit_J)
    // * }
    // */
    // while(iter != models.end()){
    //   iter2 = models.begin(); 
    //   while(iter2 != models.end()){
    //     if(iter != iter2){
    //       double m1_index = std::distance(models.begin(),iter);
    //       double m2_index = std::distance(models.begin(),iter2);
    //       
    //       if(count_params(*iter) > count_params(*iter2) && results(m1_index,0) >  results(m2_index,0)){
    //         if(results(m1_index,2) < results(m2_index,2)){
    //           results(m1_index,1) = results(m2_index,1) + fabs(R::rnorm(0, sqrt(results(m2_index,0)/10) ));
    //         }
    //       }
    //     }
    //     iter2++;
    //   }
    //   iter++; 
    // }
    
  }
  
  results.col(2) = results.col(0) + results.col(1);
  
  // Sort matrix
  arma::uvec sort_ind = sort_index(results.col(2));
  
  // Pick the "best" model
  unsigned int best_model_id = sort_ind(0);
  
  // Get the sort column index
  arma::vec ex_sort = arma::conv_to<arma::vec>::from(sort_ind);
  
  // Sort the result matrix by Criterion
  arma::mat sorted = sort_mat(results, 2);
  
  // ---------
  
  // Set up output feed
  arma::field< arma::field<arma::mat> > out(2);
  
  // Build the results matrix
  arma::field<arma::mat> ms(2);
  ms(0) = sorted;
  ms(1) = ex_sort + 1; // Sorted vector (so R can know desc IDs) ALSO change index
  
  // Create m 
  
  // If the output object is not the full model, let's inject the new results.
  if(best_model_id != full_model_index){
    
    arma::field<arma::mat> m;
    m = model_results(best_model_id);
    mod_output(0) = m(0);
    mod_output(1) = m(1);
    mod_output(5) = m(2);
    mod_output(8) = m(3);
    mod_output(9) = m(4);
    mod_output(10) = m(5);
  }
  
  
  // Output to me!
  out(0) = ms;
  out(1) = mod_output;
  
  
  return out;
}

//' @title Find the Rank Models result
//' @description Provides the core material to create an S3 object for rank.models
//' @param data A \code{vec} of data.
//' @param model_str A \code{vector<vector<string>>} that gives a list of models to test.
//' @param full_model A \code{vector<string>} that contains the largest / full model.
//' @param alpha A \code{double} that indicates the alpha level for CIs.
//' @param compute_v A \code{string} indicating the type of V matrix to generate
//' @param model_type A \code{string} that describes the model generation / transformation: 'ssm' or 'imu'
//' @param K A \code{int} that controls how many times the GMWM is run. 
//' @param H A \code{int} that controls how many bootstraps occur.
//' @param G A \code{int} that controls how many guesses occur.
//' @param robust A \code{bool} that indicates whether to use classical or robust wavelet variance.
//' @param eff A \code{double} that indicates the efficiency to use.
//' @param bs_optimism A \code{bool} that indicates whether the model selection score should be calculated with bootstrap or asymptotics.
//' @param seed A \code{unsigned int} that is the seed one wishes to use. 
//' @return A \code{field<field<field<mat>>>} that contains the model score matrix and the best GMWM model object.
//' @keywords internal
// [[Rcpp::export]]
arma::field< arma::field<arma::field<arma::mat> > >  rank_models_cpp(arma::vec& data,
                                                                     const std::vector<std::vector < std::string > >& model_str, 
                                                                     const std::vector< std::string >&  full_model,
                                                                     double alpha, 
                                                                     std::string compute_v, std::string model_type, 
                                                                     unsigned int K, unsigned int H, unsigned int G, 
                                                                     bool robust, double eff, bool bs_optimism, unsigned int seed){
  
  
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
    robust, eff, seed);
  
  return h;
}

//' @title Find the auto imu result
//' @description Provides the core material to create an S3 object for auto.imu
//' @param data A \code{mat} containing multiple columns of independent data with the same number of observations.
//' @param model_str A \code{vector<vector<string>>} that gives a list of models to test.
//' @param full_model A \code{vector<string>} that contains the largest / full model.
//' @param alpha A \code{double} that indicates the alpha level for CIs.
//' @param compute_v A \code{string} indicating the type of V matrix to generate
//' @param model_type A \code{string} that describes the model generation / transformation: 'ssm' or 'imu'
//' @param K A \code{int} that controls how many times the GMWM is run. 
//' @param H A \code{int} that controls how many bootstraps occur.
//' @param G A \code{int} that controls how many guesses occur.
//' @param robust A \code{bool} that indicates whether to use classical or robust wavelet variance.
//' @param eff A \code{double} that indicates the efficiency to use.
//' @param bs_optimism A \code{bool} that indicates whether the model selection score should be calculated with bootstrap or asymptotics.
//' @param seed A \code{unsigned int} that is the seed one wishes to use. 
//' @return A \code{field<field<field<mat>>>} that contains the model score matrix and the best GMWM model object.
//' @keywords internal
// [[Rcpp::export]]
arma::field< arma::field<arma::field<arma::mat> > >  auto_imu_cpp(arma::mat& data,
                                                                  const arma::mat& combs,
                                                                  const std::vector< std::string >&  full_model,
                                                                  double alpha, 
                                                                  std::string compute_v, std::string model_type, 
                                                                  unsigned int K, unsigned int H, unsigned int G, 
                                                                  bool robust, double eff, bool bs_optimism, unsigned int seed){
  
  
  
  
  // Number of columns to process
  unsigned int V = data.n_cols;
  
  // Create a set of unique models.
  std::set<std::vector<std::string > > models = build_model_set(combs, full_model);
  
  arma::field< arma::field<arma::field<arma::mat> > > h(V);
  
  for(unsigned int i = 0; i < V; i++){
    Rcpp::Rcout << "Generating models for the " << i + 1 << " column in the data set " << std::endl << std::endl;
    
    arma::vec signal_col = data.col(i);
    
    h(i) = model_select(signal_col,
                        models,
                        full_model,
                        model_type,
                        bs_optimism,
                        alpha,
                        compute_v, 
                        K, H, G, 
                        robust, eff, seed);
    
    Rcpp::Rcout << std::endl;
  }
  
  return h;
}
