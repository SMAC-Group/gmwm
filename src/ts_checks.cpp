#include <RcppArmadillo.h>

#include "ts_checks.h"

// Include polyroot for invertibility
#include "polyroot.h"

// Complex tools
#include "complex_tools.h"
using namespace Rcpp;


// [[Rcpp::export]]
double minroot(const arma::cx_vec& x){
  return min(
    // sqrt(x^2 + y^2)
    Mod_cpp(
      // Return roots
      do_polyroot_arma(x)
    )
  );
}

//' @title Check Invertibility Conditions
//' @description Checks the invertiveness of series of coefficients.
//' @param x A \code{cx_vec} that has a 1 appended before the coefficents. (e.g. c(1, x))
//' @return True (if outside unit circle) || False (if inside unit circle)
// [[Rcpp::export]]
bool invert_check(const arma::vec& x){
  return minroot(arma::conv_to<arma::cx_vec>::from(x)) > 1; // Outside the unit circle
}


int map_acc(int lhs, const std::pair<std::string, int> & rhs)
{
  return lhs + rhs.second;
}

int calc_map_sum(const std::map<std::string, int>& m){
  return std::accumulate(m.begin(), m.end(), 0, map_acc);  
}


//' @title Count Models
//' @description Count the amount of models that exist.
//' @param desc A \code{vector<string>} that contains the model's components.
//' @return A \code{map<string, int>} containing how frequent the model component appears.
//' @examples
//' #TBA
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

// [[Rcpp::export]]
arma::vec order_AR1s(arma::vec theta, const std::vector<std::string>& desc, const arma::field<arma::vec> objdesc){
  int AR1_old_loc = -1;
  
  double AR1_phi_prev = 0;
  
  double AR1_phi_act = 0;
  
  unsigned int i_theta = 0;
  
  for(unsigned int i = 0; i < desc.size(); i++){
    std::string element_type = desc[i];
    
    if(element_type == "AR1"){
      // Is this the first AR1 element in the stack?
      if(AR1_old_loc != -1){
        
        AR1_phi_act = theta(i_theta);
        
        // Make the largest phi value first.
        if(AR1_phi_prev < AR1_phi_act){
          
          // Put old phi in current position
          theta(i_theta) = AR1_phi_prev;
          
          // Move large phi value to old location 
          theta(AR1_old_loc) = AR1_phi_act;
          
          // Extract new sig in current position
          AR1_phi_prev = theta(i_theta+1);
          
          // Update sigma2 of the new location with the old value
          theta(i_theta + 1) = theta(AR1_old_loc+1);
    
          // Update old location with old sig2. 
          theta(AR1_old_loc+1) = AR1_phi_prev;
          
          // Store new low theta!
          AR1_phi_prev = theta(i_theta);
        }
        
        // Else: The current one is less than the previous.
        
        // Update old AR1 location 
        AR1_old_loc = i_theta;
          
      }else{ // First element, initialize values.
        AR1_old_loc = i_theta;
        AR1_phi_prev = theta(i_theta);
      }
      
      i_theta += 2;
    }else if(element_type == "ARMA"){
      i_theta += sum(objdesc(i));  
    }else{
      i_theta++;
    }
  }
  
  return theta;
}
