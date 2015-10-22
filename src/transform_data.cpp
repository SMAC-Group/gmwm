#include <RcppArmadillo.h>
using namespace Rcpp;

#include "transform_data.h"

//' @title Pseudo Logit Inverse Function
//' @description This function computes the pseudo inverse of a logit transformation of the parameters in order to constrain them to a positive domain 
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @examples
//' x.sim = rnorm(100)
//' pseudo_logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit_inv(const arma::vec& x){
  return 2.0/(1 + exp(-x)) - 1.0;
}


double pseudo_logit_inv(double x){
  return 2.0/(1 + exp(-x)) - 1.0;
}


//' @title Logit Inverse Function
//' @description This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @examples
//' x.sim = rnorm(100)
//' logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec logit_inv(const arma::vec& x){
  return 1.0/(1.0 + exp(-x));
}

double logit_inv(double x){
  return 1.0/(1.0 + exp(-x));
}

//' @title Pseudo Logit Function
//' @description This function compute the link function to constrain parameters to a positive domain.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @keywords internal
//' @examples
//' x.sim = runif(100)
//' pseudo_logit(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit(const arma::vec& x){
  arma::vec p = (x+1.0)/2.0;
  return log(p/(1.0 - p));
}

double pseudo_logit(double x){
  double p = (x+1.0)/2.0;
  return log(p/(1.0 - p));
}

//' @title Logit Function
//' @description This function computes the logit link function.
//' @param x A \code{vec} containing probabilities (e.g. -1 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @keywords internal
//' @examples
//' x.sim = runif(100)
//' logit(x.sim)
// [[Rcpp::export]]
arma::vec logit(const arma::vec& x){
  return log(x/(1.0 - x));
}

double logit(double x){
  return log(x/(1.0 - x));
}



//' @title Logit Function
//' @description This function computes the logit link function.
//' @param x A \code{vec} containing probabilities (e.g. -1 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @keywords internal
//' @examples
//' x.sim = runif(100)
//' logit(x.sim)
// [[Rcpp::export]]
arma::vec logit2(const arma::vec& x){
  double b = 2.0;
  double a = 4.0;
  
  return log((x+b)/(a - b - x));
}

double logit2(double x){
  double b = 2.0;
  double a = 4.0;
  
  return log((x+b)/(a - b - x));
}

//' @title Logit2 Inverse Function
//' @description This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @examples
//' x.sim = rnorm(100)
//' logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec logit2_inv(const arma::vec& x){
  return 4.0/(1 + exp(-x)) - 2.0;
}

double logit2_inv(double x){
  return 4.0/(1 + exp(-x)) - 2.0;
}

//' @title Transform Values for Optimization
//' @description Transform parameter guesses prior to estimating with GMWM
//' @template tsobj_cpp
//' @param model_type A \code{string} that contains the model type: \code{"imu"} or \code{"ssm"}
//' @return A \code{vec} containing the transformed guesses.
// [[Rcpp::export]]
arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
  arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
  unsigned int i_theta = 0;
  unsigned int num_desc = desc.size();
  
  
  for(unsigned int i = 0; i < num_desc; i++){
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1" ){
      
      // Apply model specific parameter transformation
      if(model_type == "imu"){
        starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
      }else{ // O.W. SSM case
        starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
      }
      
      // Increase index for phi term
      ++i_theta;
    }
    else if(element_type == "ARMA" ) {
      
      // Obtain the ARMA Model information
      arma::vec arma_params = objdesc(i);
      
      // Get the count of parameter values
      unsigned int p = arma_params(0); // AR(P)
      unsigned int q = arma_params(1); // MA(Q)
      
      
      // Determine transformation to apply
      if(model_type == "imu"){
        // All parameters but sigma
        unsigned int param_est = i_theta + p + q;
        
        // Use the pseudo logit function (0 <= x <= 1) => R
        starting.rows(i_theta, param_est - 1) = pseudo_logit(theta.rows(i_theta, param_est - 1));
        
        // Increment index
        i_theta = param_est;
      }
      else{ 
        // SSM Case
        // Two different transformations.
        if(p == 1){
          starting(i_theta) = pseudo_logit(theta(i_theta));
        }
        else if(p > 1){
          starting.rows(i_theta, i_theta + p - 1) = logit2(theta.rows(i_theta, i_theta + p - 1));
        }
        
        // Increment count by p
        i_theta += p;
        
        // Two different transformations.
        if(q == 1){
          starting(i_theta) = pseudo_logit(theta(i_theta));
        }else if(q > 1){
          starting.rows(i_theta, i_theta + q - 1) = logit2(theta.rows(i_theta, i_theta + q - 1));
        }
        
        
        // Increment count by q
        i_theta += q;
        

      }// end else
      
    }// end if
    
    if(element_type != "DR" ){
      // Take care of SIGMA2 term
      starting(i_theta) = log(theta(i_theta));
      
    }else{
      starting(i_theta) = theta(i_theta);
    } // we assume that drift in this case is identity
    
    ++i_theta;
  } // end for

  return starting;
}

//' @title Revert Transform Values for Display
//' @description Undo the previous transform of parameter guesses to obtain the GMWM estimates.
//' @template tsobj_cpp
//' @param model_type A \code{string} that contains the model type: \code{"imu"} or \code{"ssm"}
//' @return A \code{vec} containing the undone transformation of parameters.
// [[Rcpp::export]]
arma::colvec untransform_values(const arma::vec& theta, 
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
    arma::colvec result  = arma::zeros<arma::colvec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();

    for(unsigned int i = 0; i < num_desc; i++){
      
      std::string element_type = desc[i];
      // AR 1
  	  if(element_type == "AR1"){
        if(model_type == "imu"){
          result(i_theta) = arma::as_scalar(logit_inv(theta.row(i_theta)));
        }else{ // ssm
          result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
        }
  	    ++i_theta;
  	  }
      else if(element_type == "ARMA" ) {
        
        // Obtain the ARMA Model information
        arma::vec arma_params = objdesc(i);
        
        // Get the count of parameter values
        unsigned int p = arma_params(0); // AR(P)
        unsigned int q = arma_params(1); // MA(Q)
        

        // Determine transformation to apply
        if(model_type == "imu"){
          // All parameters but sigma
          unsigned int param_est = i_theta + p + q;
          
          // Use the invs pseudo logit function R => (0 <= x <= 1) 
          result.rows(i_theta, param_est - 1) = pseudo_logit_inv(theta.rows(i_theta, param_est - 1));
          
          // Increment index
          i_theta = param_est;

        }else{ 
            // SSM Case
            
            // Two different transformations.
            if(p == 1){
              result(i_theta) = pseudo_logit_inv(theta(i_theta));
            }
            else if(p > 1){
              result.rows(i_theta, i_theta + p - 1) = logit2_inv(theta.rows(i_theta, i_theta + p - 1));
            }
            
            // Increment count by p
            i_theta += p;
            
            // Two different transformations.
            if(q == 1){
              result(i_theta) = pseudo_logit_inv(theta(i_theta));
            }else if(q > 1){
              result.rows(i_theta, i_theta + q - 1) = logit2_inv(theta.rows(i_theta, i_theta + q - 1));
            }
            
            // Increment count by q
            i_theta += q;
          }// end else
          
      } // end if

      if(element_type != "DR" ){
        // Take care of SIGMA2 term or other term...
        result(i_theta) = exp(theta(i_theta));
        
      }else{
        result(i_theta) = theta(i_theta);
      } // we assume that drift in this case is identity
      
      ++i_theta;
  }  

  return result;
}