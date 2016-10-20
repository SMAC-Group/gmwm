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
#include "guess_values.h"

// Transform data using transform_values
#include "transform_data.h"

// Use yannick's starting obj
#include "objective_functions.h"

// Inline function for square
#include "inline_functions.h"

// Include sampler for randomization of draws
#include "sampler.h"

// Include invertibility check.
#include "ts_checks.h"

// Include ARMAtoMA_cpp
#include "rtoarmadillo.h"

// for estimator domination
#include "process_to_wv.h"


// ------------- New

// What estimator dominates the other? 
// see sim file for formula
// not used anymore. see
std::string dom_process(double first_wv, double ci_low, double ci_high){
  
  arma::vec s(2);
  double logof2 = log(2);
  s(0) = log(ci_high/first_wv)/logof2 /2.0;
  s(1) = log(ci_low/first_wv)/logof2 /2.0;
  
  if(arma::max(s) < -.5){
    return "QN";
  }else if(arma::min(s) > -.5){
    return "AR1";
  }
  
  return "WN";
}

// Obtain the drift slope from data


double dr_slope(const arma::vec& data){
  return (data.max()-data.min())/ double(data.n_elem);
}



// ------------- Draw functions

double draw_rw(double sigma2_total, int N){
  // sigma^2/(N*10^5), sigma^2 / N
  return R::runif(0.00001*sigma2_total/double(N), sigma2_total/double(N));
}

double draw_qn_dom(double sigma2_total){
  return R::runif(sigma2_total/8.0, sigma2_total/3.0);
}

// sigma^2 /2 * 1/(10^5), sigma^2/2 * 2/100
double draw_qn_weak(double sigma2_total){
  return R::runif(sigma2_total*0.000005, sigma2_total/100.0);
}

double draw_wn_dom(double sigma2_total){
  return R::runif(sigma2_total/2.0, sigma2_total);
}

// sigma^2/10^5
double draw_wn_weak(double sigma2_total){
  return R::runif(0.00001*sigma2_total, .1*sigma2_total);
}

double draw_drift(double ranged){
  return R::runif(ranged/100.0, ranged/2.0);
}

arma::vec draw_ar1(double sigma2_total){
  // Draw from triangle distributions for phi
  double U = R::runif(0.0, 1.0/3.0);
  
  arma::vec temp(2);
  
  // Draw for phi
  temp(0) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
  
  // 1 - phi^2
  double val = (1-square(temp(0)));
  
  // (sigma^2/2 * (1-phi^2)), (sigma^2 * (1-phi^2))
  temp(1) = R::runif(0.5*sigma2_total*val, sigma2_total*val);
  
  return temp; 
}

arma::vec draw_ar1_memory_add(double sigma2_total, double last_phi){
  arma::vec temp(2);
  
  // Draw for phi
  temp(0) = R::runif(std::max(0.95,last_phi), 0.999995);
  
  // 1 - phi^2
  double val = (1-square(temp(0)));
  
  // (sigma^2/ 10^5 * (1-phi^2)), 2*(sigma^2 * (1-phi^2))/100
  temp(1) = R::runif(0.00001*sigma2_total*val, sigma2_total*val/50.0);
  
  return temp;
}

arma::vec draw_ar1_memory(double sigma2_total, double last_phi){
  
  arma::vec temp(2);
  
  // Draw for phi
  temp(0) = R::runif(std::max(0.9,last_phi), 0.999995);
  
  // 1 - phi^2
  double val = (1-square(temp(0)));
  
  // (sigma^2/ 10^5 * (1-phi^2)), 2*(sigma^2 * (1-phi^2))/100
  temp(1) = R::runif(0.0, 0.01*sigma2_total*val );
    
 //   R::runif(0.00001*sigma2_total*val, sigma2_total*val/50.0);
  

  return temp;
}

arma::vec draw_ar1_memory_large(double sigma2_total, double last_phi){
  // Draw from triangle distributions for phi
  double Y = (1.0 - std::sqrt(1.0-3.0 * R::runif(0.0, 1.0/3.0)));
  
  arma::vec temp(2);
  
  // Draw for phi
  temp(0) = (.999995 - last_phi)*(1-Y) + last_phi;
  
  // 1 - phi^2
  double val = (1-square(temp(0)));
  
  // (sigma^2/2 * (1-phi^2)), (sigma^2 * (1-phi^2))
  temp(1) = R::runif(0.0, 0.01*sigma2_total*val );
 // R::runif(0.5*sigma2_total*val, sigma2_total*val);
  
  return temp; 
}


//' @title Randomly guess a starting parameter
//' @description Sets starting parameters for each of the given parameters. 
//' @param desc A \code{vector<string>} that contains the model's components.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @param model_type A \code{string} that indicates whether it is an SSM or sensor.
//' @param num_param An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
//' @param expect_diff A \code{double} that contains the mean of the first difference of the data
//' @param N A \code{integer} that contains the number of observations in the data.
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
//' @param double A \code{double} that contains the drift slope given by \eqn{\frac{max-min}{N}}{(Max-Min)/N}
//' @param G A \code{integer} that indicates how many random draws that should be performed.
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @keywords internal
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec guess_initial(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::mat& wv, const arma::vec& tau, double ranged, unsigned int G){
  
  // Obtain the sum of variances for sigma^2_total.
  
  arma::vec wv_empirical = wv.col(0);
  
  if(model_type=="ssm"){
    return guess_initial_old(desc, objdesc,
                   model_type, num_param, expect_diff, N,
                   wv_empirical, tau, G);
  }  
  
  double sigma2_total = arma::sum(wv_empirical);
  
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
  
  double min_obj_value = std::numeric_limits<double>::max();
  
  std::map<std::string, int> models = count_models(desc);
  
  // Check for domination between AR1, QN, and WN
  // wv(1), ci_high, ci_low
  std::string dom_value = dom_process(wv(0,0), wv(1,1), wv(1,2));
  
  // Expand the models count to set up bools
  int AR1 = models["AR1"] + models["GM"];
  int ARMA = models["ARMA"];

  bool QN = models["QN"];
  bool WN = models["WN"];
  bool RW = models["RW"];
  bool DR = models["DR"];
  
  bool only_qn = QN && !(AR1 > 0 || ARMA > 0 || RW || DR || WN);
  
  bool only_wn = WN && !(AR1 > 0 || ARMA > 0 || RW || DR || QN); 
  
  bool dom_wn = dom_value == "WN";
  
  bool dom_qn = dom_value == "QN";

  unsigned int num_desc = desc.size();
  
  // Generate B guesses of model parameters
  for(unsigned int g = 0; g < G; g++){
    unsigned int i_theta = 0;
    double last_phi = 0;
    int AR1_counter = dom_wn || (dom_qn &&  R::runif(0.0,1.0) < .75 );
    
    // Generate parameters for the model
    for(unsigned int i = 0; i < num_desc; i++){
      std::string element_type = desc[i];
      
      if(element_type == "AR1" || element_type == "GM"){
        
        AR1_counter++;
        
        // k*AR1 case
        if( AR1_counter >= 3){
          temp_theta.rows(i_theta, i_theta + 1) = draw_ar1_memory_large(sigma2_total, last_phi);
        }else if(AR1_counter == 2){
          temp_theta.rows(i_theta, i_theta + 1) = draw_ar1_memory(sigma2_total, last_phi);
        }else{ // AR1 without WN 
          temp_theta.rows(i_theta, i_theta + 1) = draw_ar1(sigma2_total);
        }
        
        last_phi = temp_theta(i_theta);
        
        // Increment i_theta position for 1 parameter (phi)
        // Second shift at end. (sigma2)
        i_theta++;
        
      } else if(element_type == "WN"){ // WN
        
        if(only_wn || dom_wn){
          temp_theta(i_theta) = draw_wn_dom(sigma2_total);
        }else{
          temp_theta(i_theta) = draw_wn_weak(sigma2_total);
        }
        
      }else if(element_type == "DR"){   
      
        temp_theta(i_theta) = draw_drift(ranged);
      
      }else if(element_type == "QN"){
        
        if(only_qn || dom_qn){
          temp_theta(i_theta) = draw_qn_dom(sigma2_total);
        }else{
          temp_theta(i_theta) = draw_qn_weak(sigma2_total);
        }
    
      }else if(element_type == "RW"){
      
        temp_theta(i_theta) = draw_rw(sigma2_total, N);
      
      }
      else {
        
        // Unpackage ARMA model parameter
        arma::vec model_params = objdesc(i);
        
        // Get position numbers (AR,MA,SIGMA2)
        unsigned int p = model_params(0), q = model_params(1);
        
        // Draw samples (need extra 1 for sigma2 return)
        temp_theta.rows(i_theta, i_theta + p + q) = arma_draws(p, q, sigma2_total);
        
        i_theta += p + q; // additional +1 added at end for sigma2

        // Add seasonal guessing. 
        if( model_params.n_elem > 3 && model_params(5) != 0){
          // Get position numbers (AR,MA,SIGMA2)
          unsigned int sp = model_params(2), sq = model_params(3);
          
          temp_theta.rows(i_theta, i_theta + sp + sq) = arma_draws(sp, sq, sigma2_total);
          
          i_theta += sp + sq; // additional +1 added at end for sigma2
        }
      }
      
      i_theta ++;
    } // end for
    
    // Find the minimum object value given drawn theta
    
    // Have to transform for objFunStarting function calculation
    arma::vec tvalues = transform_values(temp_theta, desc, objdesc, model_type);
    
    // Get objective function value
    double obj = objFunStarting(tvalues, desc, objdesc, model_type, wv_empirical, tau);
    
    // Big or small vs. current?
    if(min_obj_value > obj){
      min_obj_value = obj;
      starting_theta = temp_theta;
    } //end if
    
  } // end for
  
  return starting_theta;
}







// ------------- Old




//' @title Randomly guess starting parameters for AR1
//' @description Sets starting parameters for each of the given parameters. 
//' @param draw_id An \code{unsigned int} that contains the draw principles.
//' @param last_phi A \code{double} containing the last guessed phi value.
//' @param sigma2_total A \code{double} that contains the sum of all WVs. 
//' @param model_type A \code{string} that describes the model transformation.
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @keywords internal
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec ar1_draw(unsigned int draw_id, double last_phi, double sigma2_total, std::string model_type){
  arma::vec temp(2);
  
  
  if(draw_id == 0){
    if(model_type == "imu"){
      // Draw from triangle distributions for phi
      double U = R::runif(0.0, 1.0/3.0);
      
      // Draw for phi
      temp(0) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
      temp(1) = R::runif(0.5*sigma2_total*(1-square(temp(0))), sigma2_total);
    }
    else{ // ssm
      // Draw for phi
      temp(0) = R::runif(-0.9999999999999, 0.9999999999999);
      // Draw for sigma
      temp(1) = R::runif(0.0000000000001, sigma2_total);
    }
  }
  else{
    
    if(draw_id!=1){
      // Draw for phi on i >= 3
      temp(0) = R::runif(last_phi,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    else{
      // Draw for phi on i==1
      temp(0) = R::runif(0.7,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
    }
    
    // Draw for process variance
    temp(1) = R::runif(0.0, 0.01*sigma2_total*(1-square(temp(0))) ); // VERIFY THIS SHOULD BE PHI VALUE!!
    
  } // end if
  
  return temp;
}

//' @title Randomly guess starting parameters for ARMA
//' @description Sets starting parameters for each of the given parameters. 
//' @param p An \code{unsigned int} that contains the amount of AR parameters to generate.
//' @param q An \code{unsigned int} that contains the amount of MA parameters to generate.
//' @param sigma2_total A \code{double} that contains the sum of all WVs. 
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @keywords internal
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec arma_draws(unsigned int p, unsigned int q, double sigma2_total){
  // Loop index
  unsigned int i;
  
  // Estimated Sigma2
  double sigma2;
  
  // AR + MA + SIGMA2
  arma::vec arma(p+q+1);
  
  // Used for randomization (Needed prob = 0)
  arma::vec empty = arma::zeros<arma::vec>(0);
  
  // List of AR parameters
  arma::vec ar = arma::zeros<arma::vec>(p);
  
  // List of MA parameters
  arma::vec ma = arma::zeros<arma::vec>(q);
  
  // Storage for infinite MA
  arma::vec infMA(1000);
  
  // Draw start and end
  double start, end;

  // Begin drawing AR terms if they exist
  if(p != 0){
    arma::vec one = arma::ones<arma::vec>(1);
    
    // Generate AR values
    do{
      start = -.99999999, end = .999999999;
      for(i = 0; i < p; i++){
        // Draw point and move starting bounds up
        start = R::runif(start,end);
        
        // Assign picked point
        ar(i) = start; 
      }
      
    // Invertibility check we probably need to figure out a better guessing strategy...
    } while ( invert_check(arma::join_cols(one, -ar) )
                == false // not invertible.
              );

    // Randomize draws
    ar = rsample(ar, ar.n_elem, false, empty);
    
    // Export AR terms to ARMA
    arma.rows(0, p - 1) = ar;
  }
  
  
  // Reset start and end
  start = -.99999999, end = .999999999;
  
  for(i = 0; i < q; i++){
    start = R::runif(start,end);
    ma(i) = start;
  }
  
  if(q != 0){
    // Randomize
    ma = rsample(ma, ma.n_elem, false, empty);
    
    // Export the MA terms to ARMA
    arma.rows(p, p + q - 1) = ma;
  }
  
  // Obtain infinite MA process
  infMA = ARMAtoMA_cpp(ar, ma, 1000);
  
  // Obtain sigma2
  sigma2 = sigma2_total / (1 + arma::sum(arma::square(infMA)));
  
  // Store the value Do not need the -1.
  arma(p+q) = sigma2;
  
  return arma;
}


//' @title Randomly guess a starting parameter
//' @description Sets starting parameters for each of the given parameters. 
//' @param desc A \code{vector<string>} that contains the model's components.
//' @param objdesc A \code{field<vec>} that contains an object description (e.g. values) of the model.
//' @param model_type A \code{string} that indicates whether it is an SSM or sensor.
//' @param num_param An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
//' @param expect_diff A \code{double} that contains the mean of the first difference of the data
//' @param N A \code{integer} that contains the number of observations in the data.
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
//' @param B A \code{integer} that indicates how many random draws that should be performed.
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
//' @keywords internal
//' @examples
//' #TBA
// [[Rcpp::export]]
arma::vec guess_initial_old(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::vec& wv_empir, const arma::vec& tau, unsigned int B){
  
  // Obtain the sum of variances for sigma^2_total.
  double sigma2_total = arma::sum(wv_empir);
  
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
  
  double min_obj_value = std::numeric_limits<double>::max();
  
  std::map<std::string, int> models = count_models(desc);
  
  unsigned int num_desc = desc.size();
  
  unsigned int AR1_counter; // identifiability hack. =(
  double prev_phi; // ar1_draw needs external memory  
  
  // Generate B guesses of model parameters
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
    
    // Generate parameters for the model
    for(unsigned int i = 0; i < num_desc; i++){
      std::string element_type = desc[i];
      
      if(element_type == "AR1" || element_type == "GM"){
        temp_theta.rows(i_theta, i_theta + 1) = ar1_draw(AR1_counter, prev_phi, sigma2_total, model_type);
        prev_phi = temp_theta(i_theta);
        i_theta++; // needed to account for two parameters (e.g. phi + sigma2). Second shift at end.
        AR1_counter++;
      }
      else if(element_type == "WN"){  // WN
        temp_theta(i_theta) = R::runif(sigma2_total/2.0, sigma2_total);
      }
      else if(element_type == "DR"){   
        temp_theta(i_theta) = expect_diff;
      }
      else if(element_type == "QN"){
        temp_theta(i_theta) = R::runif(.0000001, sigma2_total);
      }
      else if(element_type == "RW"){
        temp_theta(i_theta) = R::runif(sigma2_total/double(N*1000.0), 2.0*sigma2_total/double(N));
      }
      else { // Unpackage ARMA model parameter
        arma::vec model_params = objdesc(i);
        
        // Get position numbers (AR,MA,SIGMA2)
        unsigned int p = model_params(0), q = model_params(1);
        
        // Draw samples (need extra 1 for sigma2 return)
        temp_theta.rows(i_theta, i_theta + p + q) = arma_draws(p, q, sigma2_total);
        
        i_theta += p + q; // additional +1 added at end for sigma2
        
        // Add seasonal guessing. 
        if( model_params.n_elem > 3 && model_params(5) != 0){
          // Get position numbers (AR,MA,SIGMA2)
          unsigned int sp = model_params(2), sq = model_params(3);
          
          temp_theta.rows(i_theta, i_theta + sp + sq) = arma_draws(sp, sq, sigma2_total);
          
          i_theta += sp + sq; // additional +1 added at end for sigma2
        }
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