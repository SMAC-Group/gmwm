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

#include "transform_data.h"

//' Pseudo Logit Inverse Function
//' 
//' This function computes the pseudo inverse of a logit transformation of the parameters in order to constrain them to a positive domain 
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @seealso \code{\link{pseudo_logit}}
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(1222)
//' 
//' # Simulate data
//' x.sim = runif(10, -1, 1)
//' 
//' # Transform
//' x.sim.transformed = pseudo_logit(x.sim)
//' 
//' # Untransform
//' x.sim.untransformed = pseudo_logit_inv(x.sim.transformed)
//' 
//' # Compare results
//' results = cbind(x.sim, x.sim.untransformed)
// [[Rcpp::export]]
arma::vec pseudo_logit_inv(const arma::vec& x){
  return 2.0/(1.0 + exp(-x)) - 1.0;
}


double pseudo_logit_inv(double x){
  return 2.0/(1.0 + exp(-x)) - 1.0;
}


//' Logit Inverse Function
//' 
//' This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @seealso \code{\link{logit}}
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(1412)
//' 
//' # Simulate data
//' x.sim = runif(10, -1, 1)
//' 
//' # Transform
//' x.sim.transformed = logit(x.sim)
//' 
//' # Untransform
//' x.sim.untransformed = logit_inv(x.sim.transformed)
//' 
//' # Compare results
//' results = cbind(x.sim, x.sim.untransformed)
// [[Rcpp::export]]
arma::vec logit_inv(const arma::vec& x){
  return 1.0/(1.0 + exp(-x));
}

double logit_inv(double x){
  return 1.0/(1.0 + exp(-x));
}

//' Pseudo Logit Function
//' 
//' This function compute the link function to constrain parameters to a positive domain.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @seealso \code{\link{pseudo_logit_inv}}
//' @keywords internal
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(1412)
//' 
//' # Simulate data
//' x.sim = runif(10, -1, 1)
//' 
//' # Transform
//' x.sim.transformed = pseudo_logit(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit(const arma::vec& x){
  arma::vec p = (x+1.0)/2.0;
  return log(p/(1.0 - p));
}

double pseudo_logit(double x){
  double p = (x+1.0)/2.0;
  return log(p/(1.0 - p));
}

//' Logit Function
//' 
//' This function computes the logit link function.
//' @param x A \code{vec} containing probabilities (e.g. -1 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @keywords internal
//' @seealso \code{\link{logit_inv}}
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(1142)
//' 
//' # Simulate data
//' x.sim = runif(10)
//' 
//' # Transform
//' x.sim.transformed = logit(x.sim)
// [[Rcpp::export]]
arma::vec logit(const arma::vec& x){
  return log(x/(1.0 - x));
}

double logit(double x){
  return log(x/(1.0 - x));
}



//' Logit2 Function
//' 
//' This function computes the logit2 link function.
//' @param x A \code{vec} containing probabilities (e.g. -2 <= x <= 2)
//' @return A \code{vec} containing logit terms.
//' @keywords internal
//' @seealso \code{\link{logit2_inv}}
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(1142)
//' 
//' # Simulate data
//' x.sim = runif(10, -2, 2)
//' 
//' # Transform
//' x.sim.transformed = logit2(x.sim)
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

//' Logit2 Inverse Function
//' 
//' This function computes the inverse of a logit transformation of the parameters.
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @keywords internal
//' @seealso \code{\link{logit2}}
//' @template author/jjb
//' @examples
//' # Set seed for reproducibility
//' set.seed(2234)
//' 
//' # Simulate data
//' x.sim = runif(10, -2, 2)
//' 
//' # Transform
//' x.sim.transformed = logit2(x.sim)
//' 
//' # Untransform
//' x.sim.untransformed = logit2_inv(x.sim.transformed)
//' 
//' # Compare results
//' results = cbind(x.sim, x.sim.untransformed)
// [[Rcpp::export]]
arma::vec logit2_inv(const arma::vec& x){
  return 4.0/(1.0 + exp(-x)) - 2.0;
}

double logit2_inv(double x){
  return 4.0/(1.0 + exp(-x)) - 2.0;
}

//' Transform Values for Optimization
//' 
//' Transform parameter guesses prior to estimating with GMWM
//' @template tsobj_cpp
//' @param model_type A \code{string} that contains the model type: \code{"imu"} or \code{"ssm"}
//' @return A \code{vec} containing the transformed guesses.
//' @template author/jjb
//' @keywords internal
// [[Rcpp::export]]
arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
  arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);

  // Count Number of Model Terms
  unsigned int num_desc = desc.size();
  
  // Indicator
  unsigned int i_theta = 0;
  
  // Cycle through transforms
  for(unsigned int i = 0; i < num_desc; i++){
    
    std::string element_type = desc[i];
    
    // AR 1
    if(element_type == "AR1" || element_type == "GM" || element_type == "MA1"){
      
      // Apply model specific parameter transformation
      if(model_type == "imu"){
        starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
      }else{ // O.W. SSM case
        starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
      }
      
      // Increase index for phi term
      ++i_theta;
    }else if(element_type == "ARMA11"){
      
      // Obtain the ARMA Model information
      arma::vec arma_params = objdesc(i);
      
      unsigned int p = arma_params(0);
      unsigned int q = arma_params(1);
      
      // Determine transformation to apply
      if(model_type == "imu"){
        
        // All parameters but sigma
        unsigned int param_est = i_theta + p + q;
        
        // Use the invs pseudo logit function R => (0 <= x <= 1) 
        starting.rows(i_theta, param_est - 1) = pseudo_logit(theta.rows(i_theta, param_est - 1));
        
        // Increment index
        i_theta = param_est;
        
      }else{ 
        // SSM Case
        
        // Two different transformations.
        if(p == 1){
          starting(i_theta) = pseudo_logit(theta(i_theta));
        }
        
        // Increment count by p
        i_theta += p;
        
        // Two different transformations.
        if(q == 1){
          starting(i_theta) = pseudo_logit(theta(i_theta));
        }
        
        // Increment count by q
        i_theta += q;
      }
    } else if(element_type == "SARIMA") {
      
      // Obtain the ARMA Model information
      arma::vec arma_params = objdesc(i);
      
      // Get the count of parameter values
      unsigned int p = arma_params(0); // AR(p)
      unsigned int q = arma_params(1); // MA(q)
      
      // Get the count of parameter values
      unsigned int sp = arma_params(2); // SAR(P)
      unsigned int sq = arma_params(3); // SMA(Q)
      
      unsigned int sfreq = arma_params(5); // season frequency
      
      // Determine transformation to apply
      if(model_type == "imu"){
        // All parameters but sigma
        unsigned int param_est = i_theta + p + q + sp + sq;
        
        // Use the invs pseudo logit function R => (0 <= x <= 1) 
        starting.rows(i_theta, param_est - 1) = pseudo_logit(theta.rows(i_theta, param_est - 1));
        
        // Increment index
        i_theta = param_est;
        
      }else{ 
        // SSM Case
        
        // Two different transformations.
        if(p == 1){
          starting(i_theta) = pseudo_logit(theta(i_theta));
        } else if(p > 1){
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
        
        if(sfreq > 0){
          // Two different transformations.
          if(sp == 1){
            starting(i_theta) = pseudo_logit(theta(i_theta));
          } else if(sp > 1){
            starting.rows(i_theta, i_theta + sp - 1) = logit2(theta.rows(i_theta, i_theta + sp - 1));
          }
          
          // Increment count by sp
          i_theta += sp;
          
          // Two different transformations.
          if(sq == 1){
            starting(i_theta) = pseudo_logit(theta(i_theta));
          }else if(sq > 1){
            starting.rows(i_theta, i_theta + sq - 1) = logit2(theta.rows(i_theta, i_theta + sq - 1));
          }
          
          // Increment count by sq
          i_theta += sq;
          
        }
      }// end else
      
      // end ARMA case 
    } 
    
    // Here we log scale the SIGMA2, RW, or QN terms
    // We are unable to identify whether a drift has a negative trend due to covariance matrix.
    // Hence, we've switched from identity to log(abs()) transform.
    if(element_type != "DR" ){
      // Take care of SIGMA2 term
      starting(i_theta) = log(theta(i_theta));
    }else{ // Handles DR term
      starting(i_theta) = log(fabs(theta(i_theta)));
    } 
    
    ++i_theta;
  } // end for
  
  return starting;
}

//' Revert Transform Values for Display
//' 
//' Undo the previous transform of parameter guesses to obtain the GMWM estimates.
//' @template tsobj_cpp
//' @param model_type A \code{string} that contains the model type: \code{"imu"} or \code{"ssm"}
//' @return A \code{vec} containing the undone transformation of parameters.
//' @template author/jjb
//' @keywords internal
// [[Rcpp::export]]
arma::vec untransform_values(const arma::vec& theta, 
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type){
  
  // Setup empty vector to store results. 
  arma::vec result  = arma::zeros<arma::vec>(theta.n_elem);
  
  // Parameter counter
  unsigned int i_theta = 0;
  
  // Number of Models
  unsigned int num_desc = desc.size();
  
  for(unsigned int i = 0; i < num_desc; i++){
    
    std::string element_type = desc[i];
    // AR 1
    if(element_type == "AR1" || element_type == "GM" || element_type == "MA1"){
      if(model_type == "imu"){
        result(i_theta) = arma::as_scalar(logit_inv(theta.row(i_theta)));
      }else{ // ssm
        result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
      }
      ++i_theta;
    } else if(element_type == "ARMA11"){
      
      // Obtain the ARMA Model information
      arma::vec arma_params = objdesc(i);
      
      unsigned int p = arma_params(0);
      unsigned int q = arma_params(1);
      
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
        
        // Increment count by p
        i_theta += p;
        
        // Two different transformations.
        if(q == 1){
          result(i_theta) = pseudo_logit_inv(theta(i_theta));
        }
        
        // Increment count by q
        i_theta += q;
      }
    } else if(element_type == "SARIMA") {
      
      // Obtain the ARMA Model information
      arma::vec arma_params = objdesc(i);
      
      // Get the count of parameter values
      unsigned int p = arma_params(0); // AR(p)
      unsigned int q = arma_params(1); // MA(q)
      
      // Get the count of parameter values
      unsigned int sp = arma_params(2); // SAR(P)
      unsigned int sq = arma_params(3); // SMA(Q)
      
      unsigned int sfreq = arma_params(5); // season frequency
      
      // Determine transformation to apply
      if(model_type == "imu"){
        // All parameters but sigma
        unsigned int param_est = i_theta + p + q + sp + sq;
        
        // Use the invs pseudo logit function R => (0 <= x <= 1) 
        result.rows(i_theta, param_est - 1) = pseudo_logit_inv(theta.rows(i_theta, param_est - 1));
        
        // Increment index
        i_theta = param_est;
        
      }else{ 
        // SSM Case
        
        // Two different transformations.
        if(p == 1){
          result(i_theta) = pseudo_logit_inv(theta(i_theta));
        } else if(p > 1){
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
        
        if(sfreq > 0){
          // Two different transformations.
          if(sp == 1){
            result(i_theta) = pseudo_logit_inv(theta(i_theta));
          } else if(sp > 1){
            result.rows(i_theta, i_theta + sp - 1) = logit2_inv(theta.rows(i_theta, i_theta + sp - 1));
          }
          
          // Increment count by sp
          i_theta += sp;
          
          // Two different transformations.
          if(sq == 1){
            result(i_theta) = pseudo_logit_inv(theta(i_theta));
          }else if(sq > 1){
            result.rows(i_theta, i_theta + sq - 1) = logit2_inv(theta.rows(i_theta, i_theta + sq - 1));
          }

          // Increment count by sq
          i_theta += sq;
          
        }
      }// end else
      
      // end ARMA case 
    } 
    
    // The remaining terms should be SIGMA2, RW, DR, WN, or QN
    // Thus, they are inversed with exp(theta)
    result(i_theta) = exp(theta(i_theta));
    
    ++i_theta;
  }  
  
  return result;
}