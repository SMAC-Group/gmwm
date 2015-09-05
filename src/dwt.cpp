#include <RcppArmadillo.h>
#include "dwt.h"

// We use QMF, Haar, Select filter, etc.
#include "wv_filters.h"

// We use reverse vec
#include "armadillo_manipulations.h"

using namespace Rcpp;

/* --------------------- Start DWT and MODWT Functions --------------------- */


//' @title Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation. 
//' @usage dwt_cpp(x, filter_name, nlevels, boundary)
//' @param x A \code{vector} with dimensions \eqn{N\times 1}{N x 1}. 
//' @param filter_name A \code{string} indicating the filter.
//' @param nlevels An \code{integer}, \eqn{J}, indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm
//' @author JJB
//' @examples
//' set.seed(999)
//' x = rnorm(2^8)
//' dwt_cpp(x, "haar", 4, boundary="periodic")
// [[Rcpp::export]]
arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name = "haar", 
                                   unsigned int nlevels = 4, std::string boundary = "periodic") {
                                     
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
      stop("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
    
  if(double(N)/double(tau) != floor(double(N)/double(tau))){
    stop("The supplied sample size ('x') must be divisible by 2^(nlevels). Either truncate or expand the number of samples.");
  }
  if(tau > N){
    stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");
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
  
  return y;
}



//' @title Maximum Overlap Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation
//' @usage modwt_cpp(x, filter_name, nlevels, boundary)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param filter_name A \code{string} indicating the filter.
//' @param nlevels An \code{integer} indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm.
//' Use this implementation to supply custom parameters instead of modwt(x),
//' which serves as a wrapper function.
//' @author JJB
//' @examples
//' set.seed(999)
//' x = rnorm(100)
//' modwt_cpp(x, "haar", 4, boundary="periodic")
// [[Rcpp::export]]
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name = "haar", 
                                   unsigned int nlevels = 4, std::string boundary = "periodic"){
  
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
      stop("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
  
  if(tau > N)
    stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");

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
  
  return y;
}


//' @title Removal of Boundary Wavelet Coefficients
//' @description Removes the first n wavelet coefficients.
//' @param x A \code{field<vec>} that contains the nlevel decomposition using either modwt or dwt.
//' @param wave_filter A \code{field<vec>} containing filter information. Only "haar" is implemented.
//' @param method A \code{string} to describe the mode. Choose between "modwt" and "dwt"
//' @return A \code{field<vec>} with boundary modwt or dwt taken care of.
//' @details 
//' The vectors are truncated by removing the first n wavelet coefficients. 
//' These vectors are then stored into the field that is returned.
//' Note: As a result, there are no NA's introduced and hence the na.omit is not needed.
//' @examples
//' x=rnorm(100)
//' brick_wall(modwt_cpp(x, "haar", 4, boundary="periodic"), select_filter("haar"), "modwt")
// [[Rcpp::export]]
arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  arma::field<arma::vec> wave_filter, std::string method = "modwt") 
{
    int m = as_scalar(wave_filter(0));

    for(unsigned int j = 0; j < x.n_elem; j++)
    {
        double binary_power = pow(2,j+1);

        unsigned int n = (binary_power - 1.0) * (m - 1.0);

        if (method == "dwt"){
            n = ceil((m - 2) * (1.0 - 1.0/binary_power));
        }
        arma::vec temp = x(j);
        unsigned int temp_size = temp.n_elem;
        n = std::min(n, temp_size);
        x(j) = temp.rows(n,temp_size-1);
    }
    
    return x;
}


/* -------------------------------- END DWT and MODWT Functions ---------------------- */