#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


//' @title Generate a white noise sequence
//' @description Generates a white noise sequence given sigma.
//' @param N An \code{integer} for signal length.
//' @param sigma_WN A \code{double} that contains process standard deviation.
//' @return wn A \code{vec} containing the white noise.
//' @examples
//' gen_white_noise(10, 1.5)
// [[Rcpp::export]]
arma::vec gen_white_noise(unsigned int N, double sigma_WN)
{
	arma::vec wn(N);
  
  for(unsigned int i=0; i< N; i++){
      wn(i) = R::rnorm(0.0, sigma_WN);
  }

	return wn;
}

//' @title Generate a drift
//' @description Generates a drift sequence with a given slope.
//' @param N An \code{integer} for signal length.
//' @param slope A \code{double} that contains drift slope
//' @return gd A \code{vec} containing the drift.
//' @examples
//' gen_drift(10, 8.2)
// [[Rcpp::export]]
arma::vec gen_drift(unsigned int N, double slope)
{
  arma::vec gd(N);
  gd.fill(slope);
	return cumsum(gd);
}


//' @title Generate an AR(1) sequence
//' @description Generate an AR sequence given phi and sig2.
//' @details This needs to be extended to AR(p) see \code{arima.sim} and \code{filter}.
//' @param N An \code{integer} for signal length.
//' @param phi A \code{double} that contains autocorrection.
//'	@param sig2 A \code{double} containing the residual variance.
//' @return gm A \code{vec} containing the AR(1) process.
//' @examples
//' gen_ar(10, 5, 1.2)
// [[Rcpp::export]]
arma::vec gen_ar1(unsigned int N, double phi, double sig2)
{

	arma::vec wn = gen_white_noise(N, sqrt(sig2));
	arma::vec gm = arma::zeros<arma::vec>(N);
	for(unsigned int i=1; i < N; i++ )
	{		
		gm(i-1) = phi*gm(i-1) + wn(i-1);
	}

	return gm;
}

//' @title Reverse Armadillo Vector
//' @description Reverses the order of an Armadillo Vector
//' @usage reverse_vec(x)
//' @param x A \code{column vector} of length N
//' @return x A \code{column vector} with its contents reversed.
//' @details Consider a vector x=[1,2,3,4,5], the function would output x=[5,4,3,2,1].
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
//' @return A \code{field<vec>} that contains:
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
arma::field<arma::vec> haar_filter() {
  
    arma::vec L = arma::vec("2");
    
    arma::vec g = arma::vec("0.7071067811865475 0.7071067811865475");
    
    arma::vec h = qmf(g);
    
    arma::field<arma::vec> out(3);
    
    out(0)=L;
    out(1)=h;
    out(2)=g;
    
    return out;
}

//' @title Select the Wavelet Filter
//' @description Constructs the wavelet filter to be used.
//' @usage select_filter(filter_name)
//' @param filter_name A \code{String} that must receive: \code{"haar"}.
//' @return info A \code{field<vec>} that contains:
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
arma::field<arma::vec> select_filter(String filter_name = "haar")
{
  
  arma::field<arma::vec> info(3);
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
//' @return y A \code{field<vec>} that contains:
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
arma::field<arma::vec> dwt_arma(arma::vec x, String filter_name = "haar", 
                                unsigned int nlevels = 4, String boundary = "periodic") {
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
  
  if(double(N)/double(tau) != floor(double(N)/double(tau))){
    ::Rf_error("The supplied sample size ('x') must be divisible by 2^(nlevels). Either truncate or expand the number of samples.");
  }
  if(tau > N){
    ::Rf_error("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");
  }

  arma::field<arma::vec> filter_info = select_filter(filter_name);
  
  int L = arma::as_scalar(filter_info(0));
  arma::vec h = filter_info(1); //check the pulls
  arma::vec g = filter_info(2);
  
  arma::field<arma::vec> y(J+1);
  
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
  return y;
}



//' @title Maximum Overlap Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation
//' @usage modwt_arma(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param filter A \code{string} indicating the filter.
//' @param n.levels An \code{integer} indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{field<vec>} that contains:
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
arma::field<arma::vec> modwt_arma(arma::vec x,  String filter_name = "haar", 
                                  unsigned int nlevels = 4, String boundary = "periodic") {
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

  arma::field<arma::vec> filter_info = select_filter(filter_name);
  
  int L = arma::as_scalar(filter_info(0));
  arma::vec ht = filter_info(1); 
  arma::vec gt = filter_info(2);
    
  // modwt transform
  double transform_factor = sqrt(2);
  ht /= transform_factor;
  gt /= transform_factor;

  arma::field<arma::vec> y(J+1);
  
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
  return y;
}


//' @title Absolute Value or Modulus of a Complex Number Squared.
//' @description Computes the squared value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus squared for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod(z) = \sqrt{(x^2 + y^2)}}
//' This function will return: \eqn{r^2 = Mod(z)^2 = x^2 + y^2}}
//' @examples
//' Mod_squared_arma(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_squared_arma(const arma::cx_vec& x){
   return pow(real(x),2) + pow(imag(x),2);
}

//' @title Absolute Value or Modulus of a Complex Number.
//' @description Computes the value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod(z) = \sqrt{(x^2 + y^2)}}
//' @examples
//' Mod_arma(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_arma(const arma::cx_vec& x){
   return sqrt(pow(real(x),2) + pow(imag(x),2));
}

//' @title Discrete Fourier Transformation for Autocovariance Function
//' @description Calculates the autovariance function (ACF) using Discrete Fourier Transformation.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the ACF.
//' @details 
//' This implementation is 2x as slow as Rs. 
//' Two issues: 1. memory resize and 2. unoptimized fft algorithm in arma.
//' Consider piping back into R and rewrapping the object. (Decrease of about 10 microseconds.)
//' @example
//' x=rnorm(100)
//' dft_acf(x)
// [[Rcpp::export]]
arma::vec dft_acf(arma::vec x){
    int n = x.n_elem;
    x.resize(2*n); // Resize fills value as zero.
    
    arma::cx_vec ff = arma::fft(x);
    
    arma::cx_vec iff = arma::conv_to< arma::cx_vec >::from( pow(real(ff),2) + pow(imag(ff),2) ); //expensive
    
    arma::vec out = arma::real( arma::ifft(iff) ) / n; // Divide out the n to normalize ifft
    
    //out.resize(n);    
    
    return out.rows(0,n-1);
}

//' @title Removal of Boundary Wavelet Coefficients
//' @description Removes the first n wavelet coefficients.
//' @param x A \code{field<vec>} that contains the nlevel decomposition using either modwt or dwt.
//' @param wave_filter A \code{field<vec>} containing filter information.
//' @param method A \code{string} to describe the mode.
//' @return A \code{field<vec>} with boundary modwt or dwt taken care of.
//' @details 
//' @example
//' x=rnorm(100)
//' brick_wall(modwt_arma(x))
// [[Rcpp::export]]
arma::field<arma::vec> brick_wall(arma::field<arma::vec> x, Rcpp::List wave_filter, String method = "modwt") 
{
    int m = wave_filter(0);
    unsigned int n = 0;
    unsigned int j = 0;
    unsigned int temp_size = 0;

    while(j < x.n_elem - 1)
    {
        double binary_power = pow(2,j+1);
        if (method == "dwt"){
            n = ceil((m - 2) * (1.0 - 1.0/binary_power));
        }
        else{
            n = (binary_power - 1.0) * (m - 1.0);
        }
        arma::vec temp = x(j);
        temp_size = temp.n_elem;
        n = std::min(n, temp_size);
        x(j) = temp.rows(n,temp_size-1);
        j++;
    }
    
    arma::vec temp = x(j);
    x(j) = temp.rows(n,temp_size-1);
    return x;
}



// [[Rcpp::export]]
arma::vec ci_eta3(arma::vec dims) {
    
    unsigned int num_elem = dims.n_elem;
    arma::vec eta3(num_elem);
    
    double binary_power = 2.0;
    for(unsigned int i = 0; i<num_elem;i++){
      eta3(i) = std::max(dims(i)/binary_power,1.0);
      binary_power *= 2.0;
    }

    return eta3;
}


//' @title Generate a Confidence intervval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a confidence interval
//' @param x A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param type A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the (1-p)*alpha confidence level 
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Lower Bounds}
//'  \item{Column 3}{Upper Bounds}
//' }
//' @details 
//' @example
//' x=rnorm(100)
//' wave_variance(brick_wall(modwt_arma(x), haar_filter()))
// [[Rcpp::export]]
arma::mat wave_variance(arma::field<arma::vec> x, String type = "eta3", double p = 0.025){
  
  unsigned int num_fields = x.n_elem;
  arma::vec x_ss(num_fields);
  arma::vec y(num_fields);
  arma::vec dims(num_fields);
  
  for(unsigned int i=0; i<num_fields;i++){
    arma::vec temp = x(i);
    dims(i) = temp.n_elem;
    x_ss(i) = dot(temp,temp);
    y(i) = x_ss(i)/dims(i);
  }
  
  arma::mat out(num_fields, 3);
  
  if(type == "eta3"){
      arma::vec eta3 = ci_eta3(dims);      
      arma::vec lower(num_fields);
      arma::vec upper(num_fields);
      
      for(unsigned int i =0; i<num_fields; i++){
        lower(i) = eta3(i) * y(i)/R::qchisq(1-p, eta3(i), 1, 0);
        upper(i) = eta3(i) * y(i)/R::qchisq(p, eta3(i), 1, 0);
      }
    
      out.col(0) = y;
      out.col(1) = lower;
      out.col(2) = upper;
  }else{
      ::Rf_error("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}
