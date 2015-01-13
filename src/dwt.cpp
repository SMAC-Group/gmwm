#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Reverse Armadillo Vector
//' @description Reverses the order of an Armadillo Vector
//' @usage reverse_vec(x)
//' @param x A \code{column vector} of length N
//' @return x A \code{column vector} with its contents reversed. 
//' @author JJB
//' @examples
//' x = 1:5
//' reverse_vec(x)
// [[Rcpp::export]]
arma::vec reverse_vec(arma::vec x) {
   std::reverse(x.begin(), x.end());
   return x;
}

//' @title Quadrature Mirror Filter
//' @description Calculate the series quadrature mirror filter (QMF). Requires a series of an even length.
//' @usage qmf(g, TRUE)
//' @param g A \code{vector} that contains the filter constants.
//' @param inverse A \code{bool} that indicates whether the inverse quadrature mirror filter is computed. 
//' By default, the inverse quadrature mirror is computed.
//' @return A \code{vector} that contains either the forward QMF (evalute in order) or the inverse QMF (reverse order). 
//' @details
//' @author
//' @examples
//' g = rep(1/sqrt(2))
//' qmf(g)
// [[Rcpp::export]]
arma::vec qmf(arma::vec g, bool inverse = true) {
  
  unsigned int L = g.n_elem;
  
  arma::vec rev_g = reverse_vec(g);
    
  for(unsigned int i = 0; i < L; i++){
  
    if( (i+!inverse) % 2 != 0){
      rev_g(i) = rev_g(i)*-1;
    }
    
  }
  
  return rev_g;
}

//' @title Haar filter construction
//' @description Creates the haar filter
//' @usage haar_filter()
//' @return A \code{list} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details
//' @author 
//' @examples
//' haar_filter()
// [[Rcpp::export]]
Rcpp::List haar_filter() {
  
    unsigned int L = 2;
    
    arma::vec g = arma::vec("0.7071067811865475 0.7071067811865475");
    
    arma::vec h = qmf(g);
    
    return Rcpp::List::create(Rcpp::Named("L",L),
                              Rcpp::Named("h",h),
                              Rcpp::Named("g",g));
}

//' @title Select the Wavelet Filter
//' @description Constructs the wavelet filter to be used.
//' @usage select_filter(filter_name)
//' @param filter_name A \code{String} that must receive: \code{"haar"}.
//' @return info A \code{list} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details 
//' The package is oriented toward using only the haar filter. If the package extends at a later time, then the supporting infrastructure is there.
//' @author JJB
//' @examples
//' select_filter("haar")
// [[Rcpp::export]]
Rcpp::List select_filter(String filter_name = "haar")
{
  
  Rcpp::List info(3);
  if(filter_name == "haar"){  
      info = haar_filter();
  }else{
      ::Rf_error("Wave Filter is not supported! See ?select_filter for supported types."); 
  }
  
  return info;
}



//' @title Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation
//' @usage dwt_arma(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param filter A \code{string} indicating the filter.
//' @param n.levels An \code{integer} indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{list} that contains:
//' \itemize{
//'  \item{" "}{}
//'  \item{""}{}
//'  \item{" "}{}
//' }
//' @details
//' 
//' @author JJB
//' @references 
//' @examples
//' set.seed(999)
// [[Rcpp::export]]
Rcpp::List dwt_arma(arma::vec x, 
               String filter_name = "haar", 
               unsigned int nlevels = 4, 
               String boundary = "periodic") {
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
      ::Rf_error("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
  
  if(double(N)/double(tau) != floor(double(N)/double(tau)))
    ::Rf_error("The supplied sample size ('x') must be divisible by 2^(nlevels). Either truncate or expand the number of samples."); 
  if(tau > N)
    ::Rf_error("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");

  Rcpp::List filter_info = select_filter(filter_name);
  
  int L = filter_info(0);
  arma::vec h = filter_info(1); //check the pulls
  arma::vec g = filter_info(2);
  
  Rcpp::List y(J+1);
  
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
  
  y(J) = x;
  y.attr("class") = "dwt";
  y.attr("boundary") = boundary;
  return y;
}



//' @title Maximum Overlap Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation
//' @usage modwt_arma(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param filter A \code{string} indicating the filter.
//' @param n.levels An \code{integer} indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"allan"}{The allan variance}
//'  \item{"errors"}{The error associated with the variance calculation.}
//' }
//' @details
//' 
//' @author JJB
//' @references 
//' @examples
//' set.seed(999)
// [[Rcpp::export]]
Rcpp::List modwt_arma(arma::vec x, 
               String filter_name = "haar", 
               unsigned int nlevels = 4, 
               String boundary = "periodic") {
  if(boundary == "periodic"){
    //
  }else if(boundary == "reflection"){
    unsigned int temp_N = x.n_elem;
    arma::vec rev_vec = reverse_vec(x);
    x.resize(2*temp_N);
    x.rows(temp_N, 2*temp_N-1) = rev_vec;
  }else{
      ::Rf_error("The supplied 'boundary' argument is not supported! Choose either periodic or reflection."); 
  }

  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2,J);
  
  if(tau > N)
    ::Rf_error("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");

  Rcpp::List filter_info = select_filter(filter_name);
  
  int L = filter_info(0);
  arma::vec ht = filter_info(1); // modwt transform
  arma::vec gt = filter_info(2);
  
  double transform_factor = sqrt(2);
  
  ht /= transform_factor;
  gt /= transform_factor;
  
  Rcpp::List y(J+1);
  
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
  
  y(J) = x;
  y.attr("class") = "dwt";
  y.attr("boundary") = boundary;
  return y;
}
