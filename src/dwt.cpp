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
#include "dwt.h"

// We use QMF, Haar, Select filter, etc.
#include "wv_filters.h"

// We use reverse vec
#include "armadillo_manipulations.h"

/* --------------------- Start DWT and MODWT Functions --------------------- */


//' @title Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation. 
//' @param x           A \code{vector} with dimensions \eqn{N\times 1}{N x 1}. 
//' @param filter_name A \code{string} indicating the filter.
//' @param nlevels     An \code{integer}, \eqn{J}, indicating the level of the decomposition.
//' @param boundary    A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @param brickwall   A \code{bool} indicating whether the a brick wall procedure should be applied to the coefficients.
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm
//' @author JJB
//' @keywords internal
//' @examples
//' set.seed(999)
//' x = rnorm(2^8)
//' dwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = TRUE)
// [[Rcpp::export]]
arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name, 
                                   unsigned int nlevels, std::string boundary, bool brickwall) {
                                     
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
      Rcpp::stop("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
    
  if(double(N)/double(tau) != floor(double(N)/double(tau))){
    Rcpp::stop("The supplied sample size ('x') must be divisible by 2^(nlevels). Either truncate or expand the number of samples.");
  }
  if(tau > N){
    Rcpp::stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");
  }

  arma::field<arma::vec> filter_info = select_filter(filter_name);
  
  int L = arma::as_scalar(filter_info(0));
  arma::vec h = filter_info(1); //check the pulls
  arma::vec g = filter_info(2);
  
  arma::field<arma::vec> y(J);
  
  for(unsigned int j = 1; j <= J; j++) {
    
    unsigned int M = N/pow(2,(j-1));
    unsigned int M_over_2 = double(M)/2;
    
    arma::vec Wj(M_over_2);
    arma::vec Vj(M_over_2);
    
    for(unsigned t = 0; t < M_over_2; t++) {
      
      int u = 2*t + 1;

      double Wjt = h(0)*x(u);
      double Vjt = g(0)*x(u);
      
      for(int n = 1; n < L; n++){
        u -= 1;
        if(u < 0){
          u = M - 1;
        } 
        Wjt += h(n)*x(u);
        Vjt += g(n)*x(u);
      }
      
      Wj[t] = Wjt;
      Vj[t] = Vjt;
    }
    
    y(j-1) = Wj;
    x = Vj;
  }
  
  
  // Apply brickwall
  if(brickwall){
    y = brick_wall(y, filter_info, "dwt");
  }
  
  return y;
}



//' @title Maximum Overlap Discrete Wavelet Transform
//' @description 
//' Calculation of the coefficients for the discrete wavelet transformation
//' @inheritParams dwt_cpp
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @keywords internal
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm.
//' Use this implementation to supply custom parameters instead of modwt(x),
//' which serves as a wrapper function.
//' @author JJB
//' @keywords internal
//' @examples
//' set.seed(999)
//' x = rnorm(100)
//' modwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = TRUE)
// [[Rcpp::export]]
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name, 
                                   unsigned int nlevels, std::string boundary, bool brickwall){
  
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
    Rcpp::stop("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
  
  if(tau > N) Rcpp::stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");

  arma::field<arma::vec> filter_info = select_filter(filter_name);
  
  int L = arma::as_scalar(filter_info(0));
  arma::vec ht = filter_info(1); 
  arma::vec gt = filter_info(2);
  
  // modwt transform
  double transform_factor = sqrt(2);
  ht /= transform_factor;
  gt /= transform_factor;

  arma::field<arma::vec> y(J);
  
  arma::vec Wj(N);
  arma::vec Vj(N);
  
  for(unsigned int j = 1; j <= J; j++) {
    for(unsigned t = 0; t < N; t++) {
      
      int k = t;

      double Wjt = ht(0)*x(k);
      double Vjt = gt(0)*x(k);
  
      for(int n = 1; n < L; n++){
        k -= pow(2, j-1);
        if(k < 0){
          k += N;
        } 
        Wjt += ht(n)*x(k);
        Vjt += gt(n)*x(k);
      }
      
      Wj[t] = Wjt;
      Vj[t] = Vjt;
    }
    
    y(j-1) = Wj;
    x = Vj;
  }
  
  // Apply brickwall
  if(brickwall){
   y = brick_wall(y, filter_info, "modwt");
  }
  
  return y;
}


//' @title Removal of Boundary Wavelet Coefficients
//' @description Removes the first n wavelet coefficients.
//' @param x           A \code{field<vec>} that contains the nlevel decomposition using either modwt or dwt.
//' @param wave_filter A \code{field<vec>} containing filter information. Only "haar" is implemented.
//' @param method      A \code{string} to describe the mode. Choose between "modwt" and "dwt"
//' @return A \code{field<vec>} with boundary modwt or dwt taken care of.
//' @keywords internal
//' @details 
//' The vectors are truncated by removing the first n wavelet coefficients. 
//' These vectors are then stored into the field that is returned.
//' Note: As a result, there are no NA's introduced and hence the na.omit is not needed.
//' @examples
//' x = rnorm(100)
//' me = modwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = FALSE)
//' brick_wall(me, select_filter("haar"), "modwt")
// [[Rcpp::export]]
arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  
                                  arma::field<arma::vec> wave_filter, 
                                  std::string method) 
{
    int m = as_scalar(wave_filter(0));

    for(unsigned int j = 0; j < x.n_elem; j++)
    {
        double binary_power = pow(2.0,double(j)+1.0);

        int n = (binary_power - 1.0) * (m - 1.0);

        if (method == "dwt"){
            n = ceil((m - 2.0) * (1.0 - 1.0/binary_power));
        }
        arma::vec temp = x(j);
        int temp_size = temp.n_elem - 1; // numbers are 0,...,(N-1)
        n = std::min(n, temp_size);
        
        // Addresses the case where all scales are removed.
        if(n != temp_size){        
          x(j) = temp.rows(n,temp_size);
        }else{
          x(j) = arma::zeros<arma::vec>(0);
        }
    }
    
    return x;
}


/* -------------------------------- END DWT and MODWT Functions ---------------------- */