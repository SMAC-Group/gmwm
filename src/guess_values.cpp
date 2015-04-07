#include <RcppArmadillo.h>
#include "guess_values.h"

// Transform data using transform_values
#include "transform_data.h"

// Use yannick's starting obj
#include "objective_functions.h"

// Inline function for square
#include "inline_functions.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec ar1_draw(unsigned int draw_id, double last_phi, double sigma_tot, std::string model_type){
  arma::vec temp(2);
  
  
  if(draw_id == 0){
    if(model_type == "imu"){
      // Draw from triangle distributions for phi
      double U = R::runif(0.0, 1.0/3.0);
      
      // Draw for phi
      temp(0) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
      temp(1) = R::runif(0.95*sigma_tot*(1-square(temp(0))), sigma_tot);
    }
    else{ // ssm
      // Draw for phi
      temp(0) = R::runif(-0.9999999999999, 0.9999999999999);
      // Draw for sigma
      temp(1) = R::runif(0.0000000000001, sigma_tot);
    }
  }
  else{
    
    if(draw_id!=1){
      // Draw for phi on i >= 3
      temp(0) = R::runif(last_phi,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    else{
      // Draw for phi on i==1
      temp(0) = R::runif(0.995,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    
    // Draw for process variance
    temp(1) = R::runif(0.0, 0.01*sigma_tot*(1-square(temp(0))) ); // VERIFY THIS SHOULD BE PHI VALUE!!
    
  } // end if
  
  return temp;
}

// [[Rcpp::export]]
unsigned int count_AR1s(std::vector<std::string> s) {
  unsigned int count = 0;

  for (unsigned int i = 0; i < s.size(); i++)
    if (s[i] == "AR1") count++;

  return count;
}

// [[Rcpp::export]]
std::map<std::string, int> count_models(const std::vector<std::string>& desc){    
  std::map<std::string, int> w;	
  
  // We want to see the only the following objects with these initial values
  w["AR1"]=0;
  w["ARMA"]=0;
  w["DR"]=0;		
  w["RW"]=0;		
  w["QN"]=0;		
  w["WN"]=0;		
  
  for (unsigned int i = 0; i < desc.size(); i++) {		
    ++w[desc[i]];		
  }		
  
  return w;		
} 

//' @title Randomly guess a starting parameter
//' @description Sets starting parameters for each of the given parameters. 
//' @param desc A \code{vector<string>} that contains the model's components.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @param model_type A \code{string} that indicates whether it is an SSM or IMU.
//' @param num_param An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
//' @param expect_diff A \code{double} that contains the mean of the first difference of the data
//' @param N A \code{integer} that contains the number of observations in the data.
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
//' @param B A \code{integer} that indicates how many random draws that should be performed.
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec guess_initial(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::vec& wv_empir, const arma::vec& tau, unsigned int B){
                          
  // Obtain the sum of variances for sigma^2_total.
  double sigma_tot = arma::sum(wv_empir);
    
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
    
  double min_obj_value = std::numeric_limits<double>::max();
  
  std::map<std::string, int> models = count_models(desc);
  
  unsigned int num_desc = desc.size();
  
  unsigned int AR1_counter; // identifiability hack. =(
  double prev_phi; // ar1_draw needs external memory  
  
  // Generate parameters for the model
  for(unsigned int b = 0; b < B; b++){
    
    unsigned int i_theta = 0;

    if(models["WN"] >= 1 && model_type=="imu"){
      AR1_counter = 2;
      prev_phi = .9;
    }
    else{ // Multiple WN in model
      AR1_counter = 0;
      prev_phi = 0;
    }
    
        
    for(unsigned int i = 0; i < num_desc; i++){
      std::string element_type = desc[i];
      
      if(element_type == "AR1"){
        temp_theta.rows(i_theta, i_theta + 1) = ar1_draw(AR1_counter, prev_phi, sigma_tot, model_type);
        prev_phi = temp_theta(i_theta);
        i_theta++; // needed to account for two parameters (e.g. phi + sigma2). Second shift at end.
        AR1_counter++;
      }
      else if(element_type == "ARMA"){
        //  This needs to be implemented.
      }
      else if(element_type == "DR"){   
        temp_theta(i_theta) = expect_diff;
      }
      else if(element_type == "QN"){
        temp_theta(i_theta) = R::runif(.0000001, sigma_tot);
      }
      else if(element_type == "RW"){
        temp_theta(i_theta) = R::runif(sigma_tot/double(N*1000.0), 2.0*sigma_tot/double(N));
      }
      else{ // WN
        temp_theta(i_theta) = R::runif(sigma_tot/2.0, sigma_tot);
      }
      i_theta ++;
    } // end for
    
    arma::vec tvalues = transform_values(temp_theta, desc, objdesc, model_type);
        
    double obj = objFunStarting(tvalues, desc, objdesc, model_type, wv_empir, tau);
    
    if(min_obj_value > obj){
      min_obj_value = obj;
      starting_theta = temp_theta;
    } //end if
  } // end for
  
  return starting_theta;
}
