#include <RcppArmadillo.h>
using namespace Rcpp;

#include "transform_data.h"

//' @title Pseudo Logit Inverse Function
//' @description This function computes the pseudo inverse of a logit transformation of the parameters in order to constrain them to a positive domain 
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @examples
//' x.sim = rnorm(100)
//' pseudo_logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit_inv(const arma::vec& x){
  return 2*exp(x)/(1 + exp(x)) -1;
}

//' @title Logit Inverse Function
//' @description This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @examples
//' x.sim = rnorm(100)
//' logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec logit_inv(const arma::vec& x){
  return 1/(1 + exp(-x));
}

//' @title Pseudo Logit Function
//' @description This function compute the link function to constrain parameters to a positive domain.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @examples
//' x.sim = runif(100)
//' pseudo_logit(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit(const arma::vec& x){
  arma::vec p = (x+1)/2;
  return log(p/(1 - p));
}

double pseudo_logit(double x){
  double p = (x+1)/2;
  return log(p/(1 - p));
}

//' @title Logit Function
//' @description This function computes the logit link function.
//' @param x A \code{vec} containing probabilities (e.g. -1 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @examples
//' x.sim = runif(100)
//' logit(x.sim)
// [[Rcpp::export]]
arma::vec logit(const arma::vec& x){
  return log(x/(1 - x));
}

double logit(double x){
  return log(x/(1 - x));
}

//' @title Transform Values for Optimization
//' @description Transform parameter guesses prior to estimating with GMWM
//' @return A \code{vec} containing the transformed guesses.
// [[Rcpp::export]]
arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
    arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // AR 1
      if(desc[i] == "AR1" ){
        
        // Apply model specific parameter transformation
        if(model_type == "imu"){
  	      starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
        }else{ // O.W. SSM case
          starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
        }
        
        // Increase index for phi term
  	    ++i_theta;
        
        // Handle the SIGMA2 term
  	    starting(i_theta) = log(theta(i_theta));
        
  	  }
      else if( desc[i] == "ARMA" ) {
        
        arma::vec arma_params = objdesc(i);
        
        unsigned int p = arma_params(0); // AR(P)
        unsigned int q = arma_params(1); // MA(Q)
        
        unsigned int param_space = i_theta + p + q;
        
        // Change pseudo_logit to logit based on model type??
        starting.rows(i_theta, param_space - 1) = pseudo_logit(theta.rows(i_theta, param_space - 1));
        
        // Increment index
        i_theta += param_space;
        
        // Take care of SIGMA2 term
        starting(param_space) = log(theta(param_space));

      }
  	  else{
        starting(i_theta) = log(theta(i_theta));
  	  }
      ++i_theta;
  }  
    
  return starting;
}

//' @title Revert Transform Values for Display
//' @description Undo the previous transform of parameter guesses to obtain the GMWM estimates.
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
  	    result(i_theta) = exp(theta(i_theta));
  	  }
      else if(element_type == "ARMA" ) {
        
        arma::vec arma_params = objdesc(i);
        
        unsigned int p = arma_params(0); // AR(P)
        unsigned int q = arma_params(1); // MA(Q)
        
        unsigned int param_space = i_theta + p + q;
        
        // Change pseudo_logit_inv to logit_inv based on model type??
        result.rows(i_theta, param_space - 1) = pseudo_logit_inv(theta.rows(i_theta, param_space - 1));
        
        // Increment index to account for P+Q values
        i_theta += param_space;
        
        // Take care of SIGMA2 term
        result(param_space) = exp(theta(param_space));

      }
  	  else {
        result(i_theta) = exp(theta(i_theta));
  	  }
      
      ++i_theta;
  }  
    
  return result;
}
