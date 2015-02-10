#include <RcppArmadillo.h>
#include <string>
#include <map>
// #include <omp.h>
using namespace arma;
using namespace Rcpp;

// // [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::vec diff_arma(const arma::vec& x, unsigned int lag = 1){
  unsigned int n=x.n_elem;
  return (x.rows(lag,n-1) - x.rows(0,n-lag-1));
}

//' @title Generate a white noise sequence
//' @description Generates a white noise sequence given sigma.
//' @param N An \code{integer} for signal length.
//' @param sigma_WN A \code{double} that contains process standard deviation.
//' @return wn A \code{vec} containing the white noise.
//' @examples
//' gen_wn(10, 1.5)
// [[Rcpp::export]]
arma::vec gen_wn(const unsigned int N, const double sigma_WN)
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
//' gen_dr(10, 8.2)
// [[Rcpp::export]]
arma::vec gen_dr(const unsigned int N, const double slope)
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
arma::vec gen_ar1(const unsigned int N, const double phi, const double sig2)
{

	arma::vec wn = gen_wn(N, sqrt(sig2));
	arma::vec gm = arma::zeros<arma::vec>(N);
	for(unsigned int i=1; i < N; i++ )
	{		
		gm(i-1) = phi*gm(i-1) + wn(i-1);
	}

	return gm;
}

//' @title Generate a random walk without drift
//' @description Generates a random walk without drift.
//' @param N An \code{integer} for signal length.
//' @return grw A \code{vec} containing the random walk without drift.
//' @examples
//' gen_rw(10, 8.2)
// [[Rcpp::export]]
arma::vec gen_rw(const unsigned int N)
{
  arma::vec grw(N);
  grw.imbue( norm_rand );
  return cumsum(grw);
}



//' @title Logit Inverse Function
//' @description This function computes the probabilities
//' @param x A \code{vec} containing real numbers.
//' @return A \code{vec} containing logit probabilities.
//' @example
//' x.sim = rnorm(100)
//' pseudo_logit_inv(x.sim)
// [[Rcpp::export]]
arma::vec pseudo_logit_inv(const arma::vec& x){
  return 2*exp(x)/(1 + exp(x)) -1;
}

//' @title Logit Function
//' @description This function compute the link term.
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
//' @return A \code{vec} containing logit terms.
//' @example
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


//' @title White Noise to WV
//' @description This function compute the WV (haar) of a White Noise process
//' @param sig2 A \code{double} corresponding to variance of WN
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the white noise.
//' @example
//' x.sim = cumsum(rnorm(100000))
//' waveletVariance(x.sim)
//' wv.theo = wn_to_wv(1, floor(log(length(x.sim),2)))
//' lines(wv.theo$scales,wv.theo$WV, col = "red")
// [[Rcpp::export]]
arma::vec wn_to_wv(double sig2, arma::vec Tau){
	return sig2/Tau;
}


//' @title Random Walk to WV
//' @description This function compute the WV (haar) of a Random Walk process
//' @param sig2 A \code{double} corresponding to variance of RW
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the random walk.
//' @example
//' x.sim = cumsum(rnorm(100000))
//' waveletVariance(x.sim)
//' wv.theo = rw_to_wv(1,x.sim)
//' lines(wv.theo$scales,wv.theo$WV, col = "red")
// [[Rcpp::export]]
arma::vec rw_to_wv(double sig2, const arma::vec& Tau){
	return sig2*((2*arma::square(Tau) + 2)/(24*Tau));
}


//' @title Drift to WV
//' @description This function compute the WV (haar) of a Drift process
//' @param omega A \code{double} corresponding to variance of drift
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the drift.
//' @example
//' x.sim = 1:1000
//' waveletVariance(x.sim)
//' wv.theo = dr_to_wv(1,floor(log(length(x.sim),2)))
//' lines(wv.theo$scales,wv.theo$WV, col = "red")
// [[Rcpp::export]]
arma::vec dr_to_wv(double omega,const arma::vec& Tau){
	return (omega*omega)*arma::square(Tau)/16;
}

//' @title AR1 process to WV
//' @description This function compute the WV (haar) of an AR(1) process
//' @param omega A \code{double} corresponding to variance of drift
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the drift.
//' @example
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sig2 = 4 )
//' waveletVariance(x.sim)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, floor(log(length(x.sim),2)))
//' lines(wv.theo$scales,wv.theo$WV, col = "red")
// [[Rcpp::export]]
arma::vec ar1_to_wv(double phi, double sig2, const arma::vec& Tau){
  unsigned int size_tau = Tau.n_elem;
  arma::vec temp_term(size_tau);
  arma::vec temp_term_redux(size_tau);
  for(unsigned int i=0; i< size_tau; i++){
    temp_term(i) = 4*pow(phi,(Tau(i)/2 + 1));
    temp_term_redux(i) = pow(phi,(Tau(i)+1));
  }
	return ((Tau/2 - 3*phi - Tau/2*pow(phi,2) + temp_term - temp_term_redux)/(arma::square(Tau/2)*pow(1-phi,2)*(1-pow(phi,2)))*sig2)/2;
}

//' @title Expected value DR
//' @description This function computes the expected value of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the expected value of the drift.
//' @example
//' # Add at a later time
// [[Rcpp::export]]
arma::vec e_drift(double omega, int n_ts){
  arma::vec out(1);
  out(0) = omega*(n_ts + 1.0)/2.0;
  return out;
}

//' @title Second moment DR
//' @description This function computes the second moment of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the second moment of the drift.
//' @example
//' # Add at a later time
// [[Rcpp::export]]
arma::vec m2_drift(double omega, int n_ts){
  arma::vec out(1);
  out(0)=(omega*omega)*(double(n_ts*n_ts)/3.0 + double(n_ts)/2.0 + 1.0/6.0);
  return out;
}

//' @title Variance DR
//' @description This function computes the variance of a drift process.
//' @param omega A \code{double} corresponding to variance of drift.
//' @param n_ts An \code{int} indicating the length of the time series.
//' @return A \code{vec} containing the variance of the drift.
//' @example
//' # Add at a later time
// [[Rcpp::export]]
arma::vec var_drift(double omega, int n_ts){
  // Compute m1
	arma::vec m1 = e_drift(omega, n_ts);
	
	// Compute m2
	arma::vec m2 = m2_drift(omega, n_ts);
	
	// Compute var
  return (m2 - m1*m1)*double(n_ts)/double(n_ts-1.0);
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
  
    arma::vec L(1);
    L(0) = 2.0;
    
    arma::vec g(2);
    g.fill(0.7071067811865475);
    
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
      stop("Wave Filter is not supported! See ?select_filter for supported types."); 
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
arma::field<arma::vec> modwt_arma(arma::vec x, String filter_name = "haar", 
                                   unsigned int nlevels = 4, String boundary = "periodic") {
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
arma::vec Mod_squared_arma(arma::cx_vec x){
   return arma::square(arma::real(x)) + arma::square(arma::imag(x));
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
arma::vec Mod_arma(arma::cx_vec x){
   return arma::sqrt(arma::square(arma::real(x)) + arma::square(arma::imag(x)));
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
    
    arma::cx_vec iff = arma::conv_to< arma::cx_vec >::from( arma::square(arma::real(x)) + arma::square(arma::imag(x)) ); //expensive
    
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
arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  arma::field<arma::vec> wave_filter, String method = "modwt") 
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


//' @title Generate eta3 confidence interval
//' @description Computes the eta3 CI
//' @param y A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param dims A \code{String} indicating the confidence interval being calculated.
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
arma::mat ci_eta3(arma::vec y,  arma::vec dims, double p) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);

    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2,i+1),1.0);
      out(i,1) = eta3 * y(i)/R::qchisq(1-p, eta3, 1, 0); // Lower CI
      out(i,2) = eta3 * y(i)/R::qchisq(p, eta3, 1, 0); // Upper CI
    }

    out.col(0) = y;

    return out;
}


//' @title Generate a Confidence intervval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a chi-squared confidence interval
//' @param x A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param type A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the (1-p)*alpha confidence level 
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @details 
//' @example
//' x=rnorm(100)
//' wave_variance(brick_wall(modwt_arma(x), haar_filter()))
// [[Rcpp::export]]
arma::mat wave_variance( arma::field<arma::vec> x, String type = "eta3", double p = 0.025){
  
  unsigned int num_fields = x.n_elem;
  arma::vec y(num_fields);
  arma::vec dims(num_fields);
  
  for(unsigned int i=0; i<num_fields;i++){
    arma::vec temp = x(i);
    dims(i) = temp.n_elem;
    y(i) = dot(temp,temp)/dims(i);
  }
  
  arma::mat out(num_fields, 3);
  
  if(type == "eta3"){
      out = ci_eta3(y,dims,p);      
  }
  else if(type == "none"){
      out.col(0) = y;
  }
  else{
      stop("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}

//' @title Computes the (MODWT) wavelet variance
//' @description Calculates the (MODWT) wavelet variance
//' @param x A \code{vec} that contains the signal
//' @param strWavelet A \code{String} indicating the type of wave filter to be applied. Must be "haar"
//' @param compute_v A \code{String} that indicates covariance matrix multiplication. 
//' @return A \code{field<mat>} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance},
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//'   \item{"wavelet"}{Filter Used}
//'   \item{"scales"}{Scales}
//'   \item{"V"}{Asymptotic Covariance Matrix}
//'   \item{"up_gauss"}{Upper Gaussian CI}
//'   \item{"dw_gauss"}{Lower Guassian CI}
//' }
//' @details 
//' The underlying code should be rewritten as a class for proper export.
//' @example
//' x=rnorm(100)
//' wavelet_variance_arma(x, "haar", "diag")
// [[Rcpp::export]]
arma::field<arma::mat> wavelet_variance_arma(const arma::vec& signal, String strWavelet="haar", String compute_v = "no") {

  // Set p-value for (1-p)*100 ci
  double p = 0.025;
  
  // Length of the time series
  unsigned int n_ts = signal.n_elem;
  
  // Compute number of scales considered
  unsigned int nb_level = floor(log2(n_ts));
  
  // MODWT transform
  arma::field<arma::vec> signal_modwt = modwt_arma(signal, strWavelet, nb_level);
  arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, select_filter(strWavelet));
  
  // Compute wavelet variance  
  arma::mat vmod = wave_variance(signal_modwt_bw, "eta3", p);
  
  // Define scales
  arma::vec scales(nb_level);
  for(unsigned int i=0; i< nb_level;i++){
    scales(i) = pow(2,i+1);
  }

  // Compute asymptotic covariance matrix
  arma::mat V = diagmat(arma::ones<arma::vec>(nb_level));
  //arma::vec up_gauss(nb_level);
  //arma::vec dw_gauss(nb_level);
  if (compute_v == "full" || compute_v == "diag"){
    if (compute_v == "full"){
      //V = compute_full_V(signal_modwt) // missing in action. Roberto's code I think had it.
    } 
    if (compute_v == "diag"){
      
      unsigned int num_field = signal_modwt.n_elem;
      
      arma::vec Aj(signal_modwt.n_elem);
      
      for(unsigned int i = 0; i < num_field; i++){
        // Autocovariance using the Discrete Fourier Transform
        arma::vec temp = dft_acf(signal_modwt(i));
        
        // Sum(V*V) - first_element^2 /2
        Aj(i) = dot(temp,temp) - temp(0)*temp(0)/2;
      }
      // Create diagnoal matrix (2 * Aj / length(modwt_d1)). Note: All modwt lengths are the same. Update if dwt is used.      
      V = diagmat(2 * Aj / signal_modwt(0).n_elem);
    }
    
    // Compute confidence intervals
    //up_gauss = vmod.col(0) + R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
    //dw_gauss = vmod.col(0) - R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
  }
  else{
    arma::vec temp = arma::ones<arma::vec>(nb_level);
    V = diagmat(temp);
    //up_gauss.fill(datum::nan);
    //dw_gauss.fill(datum::nan);
  }
  
  // export
  arma::field<arma::mat> out(5);
  out(0) = vmod.col(0); // wave variance
  out(1) = vmod.col(1); // low ci
  out(2) = vmod.col(2); // high ci
  out(3) = scales; // scales
  out(4) = V; //V
  //out(5) = up_gauss; // up_guass ci
  //out(6) = dw_gauss; // low_guass ci
  
  // Define structure "wav.var"
  return out;
}

// 
inline arma::vec theoretical_wv(const arma::vec& theta, const std::vector<std::string>& desc,
                                const arma::vec& wv_empir,
                                const arma::vec& tau, int N){
  
  unsigned int num_desc = desc.size();
  unsigned int i_theta = 0;
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
  
  for(unsigned int i = 0; i < num_desc; i++){
    // AR 1
    if(desc[i] == "AR1"){
      double phi = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
      ++i_theta;
      double sig2 = exp(theta(i_theta));
      ++i_theta;
      
      // Compute theoretical WV
      wv_theo += ar1_to_wv(phi,sig2,tau);
    }
    // RW
    else if(desc[i] == "RW"){
      double sig2 = exp(theta(i_theta));
      ++i_theta;
      wv_theo += rw_to_wv(sig2,tau);
    }
    // DR
    else if(desc[i] == "DR"){
      double drift = exp(theta(i_theta));
      ++i_theta;
      wv_theo += dr_to_wv(drift,tau);
    }
    // WN
    else{
      double sig2 = exp(theta(i_theta));
      ++i_theta;
      wv_theo += wn_to_wv(sig2,tau);
    }
    
  }

  return wv_theo;
}

// hiding this function for the moment
double objFunStarting(const arma::vec& theta, const std::vector<std::string>& desc, 
                         const arma::vec& wv_empir, const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, wv_empir, tau, N);
  arma::vec standardized = 1-wv_theo/wv_empir;
	// Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

double objFun(const arma::vec& theta, const arma::mat& omega,
                 const std::vector<std::string>& desc, const arma::vec& wv_empir,
                 const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, wv_empir, tau, N);

  // Compute quandratic form
	arma::vec dif = wv_theo - wv_empir;
	return arma::as_scalar(trans(dif)*omega*dif);
}

// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta, const std::vector<std::string>& desc,
                          const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        _["desc"] = desc,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}



// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, const std::vector<std::string>& desc, const arma::mat& omega, const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun),
                        _["omega"] = omega,
                        _["desc"] = desc,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

// [[Rcpp::export]]
arma::vec gen_model(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc){
    arma::vec x  = arma::zeros<arma::vec>(N);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
  	for(unsigned int i = 0; i < num_desc; i++){
  	  // AR 1
  	  if(desc[i] == "AR1" || desc[i] == "AR1P"){
  	    double phi = theta(i_theta);
  	    ++i_theta;
  	    double sig2 = theta(i_theta);
  	    ++i_theta;
  	    
  	    // Compute theoretical WV
  	    x += gen_ar1(N, phi, sig2);
  	  }
  	  // RW
  	  else if(desc[i] == "RW"){
  	    ++i_theta;
  	    x += gen_rw(N);
  	  }
  	  // DR
  	  else if(desc[i] == "DR"){
  	    double drift = theta(i_theta);
  	    ++i_theta;
  	    x += gen_dr(N, drift);
  	  }
  	  // WN
  	  else if(desc[i] == "WN"){
  	    double sigmaWN = theta(i_theta);
  	    ++i_theta;
  	    x += gen_wn(N, sigmaWN);
  	  }
  	  // EX
  	  else{
  	    ++i_theta;
  	  }
  }  
    
  return x;
}


// [[Rcpp::export]]
arma::vec set_starting_values(const arma::vec& theta, const std::vector<std::string>& desc){
    arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
  	  // AR 1
  	  if(desc[i] == "AR1" ){
  	    starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
  	    ++i_theta;
  	    starting(i_theta) = log(theta(i_theta));
  	    ++i_theta;
  	  }
  	  else{
        starting(i_theta) = log(theta(i_theta));
  	    ++i_theta;
  	  }
  }  
    
  return starting;
}


// [[Rcpp::export]]
arma::rowvec set_result_values(const arma::vec& theta, const std::vector<std::string>& desc){
    arma::rowvec result  = arma::zeros<arma::rowvec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // AR 1
  	  if(desc[i] == "AR1" || desc[i] == "AR1P"){
  	    result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
  	    ++i_theta;
  	    result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
  	  else if(desc[i] != "EX"){
        result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
      else{
        // EX
        ++i_theta;
  	  }
  }  
    
  return result;
}

// [[Rcpp::export]]
arma::vec gmwm_bootstrapper(const arma::vec&  theta, const std::vector<std::string>& desc, 
                            unsigned int tau, unsigned int N, 
                            unsigned int B = 100, bool var_or_mu = false){
  unsigned int nb_level = floor(log2(N));
  	
	arma::mat res(B, tau+1);
	for(unsigned int i=0; i<B; i++){
		arma::vec x = gen_model(N, theta, desc);

    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_arma(x, "haar", nb_level);
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter());
  
		arma::vec wv_x = wave_variance(signal_modwt_bw, "none").col(0);
    //waveletVariance(x, compute.v = "diag", verbose = FALSE)$variance
    
    arma::vec temp(1);
    if(var_or_mu){
      temp(0) = arma::var(x);
    }
    else{
      temp(0) = arma::mean(x);
    }
	  res.row(i) = arma::trans(join_cols(temp,wv_x));
	}
	return cov(res);
}

inline arma::vec ar1_draw(unsigned int num_ars, double sigma_tot){
  
  unsigned int num_params = 2*num_ars;
  arma::vec temp(num_params);
  
  for(unsigned int i = 0; i < num_ars; i++){
    // Draw from triangle distributions for phi
    double U = R::runif(0.0, 2.0/3.0);
    Rcpp::Rcout << "Random Draw: " << U << std::endl;
    if(U < 1.0/3.0){
      temp(2*i) = 1.0/3 - sqrt(1.0/9.0 - U/3.0);
      Rcpp::Rcout << "Case 1: " << temp(2*i) << std::endl;
    }else{
      temp(2*i) = 0.5 + sqrt(2.0/3.0*(U-1.0/3.0));
      Rcpp::Rcout << "Case 2: " << temp(2*i) << std::endl;
    }
    // Draw from uniform for 
    temp(2*i+1) = R::runif(sigma_tot/100000.0, sigma_tot);
    Rcpp::Rcout << "Sampled sigma: " << temp(2*i+1) << std::endl;
  }
  
  Rcpp::Rcout << "Temp Vector looks like: " << temp << std::endl;
      
  return temp;
}

inline arma::vec unif_sigma_sample(unsigned int num, double sigma_tot){
  arma::vec temp(num);
  
  for(int i = 0; i<num; i++){
    temp(i) = R::runif(0.0,sigma_tot);
  }
  
  return temp;
}

//' @title Randomly guess a starting parameter
//' @description Sets starting parameters for each of the given parameters. 
//' @param signal A \code{vec} that contains the data
//' @param w A \code{std::map<std::string,int>} that lists supported models and the amount in the model.
//' @param num_params An \code{unsigned int} 
//' @param compute_v A \code{String} that indicates covariance matrix multiplication. 
//' @return A \code{vec} containing smart parameter starting guesses to be iterated over.
arma::vec guess_initial(arma::vec signal, std::map< std::string ,int>& w, const std::vector<std::string>& desc,
                        unsigned int num_param, const arma::vec& wv_empir, const arma::vec& tau, 
                        unsigned int N, unsigned int B=1000){
  
    // Obtain the sum of variances for sigma^2_total.
  double sigma_tot = arma::sum(wv_empir);
  
  Rcpp::Rcout << "Sum of wv_empir for sigma^2_total: " << sigma_tot << std::endl;
  
  arma::vec starting_theta = arma::zeros<arma::vec>(num_param);
  arma::vec temp_theta = arma::zeros<arma::vec>(num_param);
    
  double min_obj_value = std::numeric_limits<double>::max();
  
  // Generate parameters for the model
  for(unsigned int b = 0; b < B; b++){
    
    unsigned int i_theta = 0;
    
    for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
      int out = p->second;
      if(out > 0){
        std::string type = p->first;
        if(type == "AR1"){
          temp_theta.rows(i_theta, i_theta + 2*out - 1) = ar1_draw(out, sigma_tot);
          i_theta += 2*out;
        }
        else if(type == "DR"){   
          temp_theta.rows(i_theta, i_theta + out - 1).fill( mean(diff_arma(signal,1))  );
          i_theta += out;
        }
        else if(type == "QN"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot);
          i_theta += out;
        }
        else if(type == "RW"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot)/signal.n_elem;
          i_theta += out;
        }
        else{ // WN
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot);
          i_theta += out;
        }
      }
    }
    Rcpp::Rcout << "Temp Theta looks like: " << temp_theta << std::endl;

    double obj = objFunStarting(temp_theta, desc, wv_empir, tau, N);
    
    Rcpp::Rcout << "Objective value is:" << obj << std::endl;
    if(min_obj_value > obj){
      Rcpp::Rcout << "Objective value " << obj << " was lower than " << min_obj_value << std::endl;
      min_obj_value = obj;
      starting_theta = temp_theta;
    }
  }
  
  Rcpp::Rcout << "Final Objective value " << min_obj_value << std::endl;

  return starting_theta;
}


//' @title Count Number of Model instances 
//' @description Return a model count
//' @param desc A \code{vec} that contains the model type
//' @return A \code{vec} with the model counts.
//' @details Hard coded value for the number of models we currently support. 
//' The types of models supported are:
//' \itemize{
//'   \item{0}{AR1}
//'   \item{1}{Drift}
//'   \item{2}{White Noise}
//'   \item{3}{Random Walk}
//'   \item{4}{Q Noise}
//' }
//' @details
//' old counting method
//' @example
//' count_models(1,2,3,4,4,0,0,0,3,4)
// [[Rcpp::export]]
arma::vec count_models_numeric(arma::vec desc) {
    arma::vec models = arma::zeros<arma::vec>(5); // Change if we add / subtract more terms 
    for (int i = 0; i < desc.n_elem; i++) {
        models(desc(i)) += 1;
    }
    return models;
}


// Counts 
inline std::map<std::string, int> counted_map(const std::vector<std::string>& desc){
  std::map<std::string, int> w;
  w["AR1"]=0;
  w["DR"]=0;
  w["RW"]=0;
  w["QN"]=0;
  w["WN"]=0;

  for (int i = 0; i < desc.size(); i++) {
        ++w[desc[i]];
  }
  
  return w;
} 

// 
inline unsigned int count_params(std::map<std::string, int>& w) {
  unsigned int num_params = 0;     
  for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
    std::string type = p->first;
    int num_models = p->second;
    if(type != "AR1" && num_models != 0){
      ++num_params;
    }
    else{
      if(num_models > 0){
        num_params += (2*num_models);
      }
    }
  }
  return num_params;
}

//' @title Count Number of Models (Alphanumeric)
//' @description Return a model count
//' @param desc A \code{vector<string>} that contains the model type
//' @return A \code{vec} with the model counts.
//' @details 
//' The types of models supported are:
//' \itemize{
//'   \item{AR1}
//'   \item{DR}
//'   \item{QN}
//'   \item{RW}
//'   \item{WN}
//' }
//' Note: This order is due to the fact that the models are created by inserting into a C++ Map. The map sorts these values...
//' @example
//' count_models_alpha(c("AR1","AR1","DR","WN"))
arma::vec count_models_alpha(const std::vector<std::string>& desc) {
  
    std::map<std::string, int> w = counted_map(desc);

    arma::vec num_models = arma::zeros<arma::vec>(w.size());
    
    for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
              int out = p->second;
              num_models(std::distance(w.begin(), p)) = out;
    }
    return num_models;
}

//' @title GMWM Step
//' 
//' @example
//' x=rnorm(100)
//' GMWM_adv
// [[Rcpp::export]]
arma::mat GMWM(const std::vector<std::string>& desc, const arma::vec& signal,
               const arma::mat& V, const arma::vec& wv_empir,
               const arma::vec& tau, unsigned int N, unsigned int B = 1000){
  
  // Count the number of models we are working with
  std::map<std::string, int> w = counted_map(desc);
  
  // Return the total number of parameters we need to setup.
  unsigned int num_param = count_params(w);
  
  Rcpp::Rcout << "Number params in model: " << num_param << std::endl;


  // Initialize it
  arma::mat GMWM(1,num_param);
  
  // Give it a guess
  arma::vec guess_me = guess_initial(signal, w, desc, num_param, wv_empir, tau, N, B);

  Rcpp::Rcout << "Initial Guess (untransformed) is: " << guess_me << std::endl;

  // Starting values
  arma::vec starting_theta = set_starting_values(guess_me, desc);
  
  Rcpp::Rcout << "Initial Guess (transformed) is: " << starting_theta << std::endl;
  
  // Optimize using the idea Yannick proposed in GMWM for INF Design.
  starting_theta = Rcpp_OptimStart(starting_theta, desc, tau, wv_empir, N);

  Rcpp::Rcout << "Optimized starting values (transformed) by Yannick are: " << starting_theta << std::endl;
      
  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
    	
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
    
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, omega, tau, wv_empir, N);
  
  Rcpp::Rcout << "GMWM estimator transformed" << estim_GMWM << std::endl;

  // Save results
  GMWM.row(0) = set_result_values(estim_GMWM, desc);

  Rcpp::Rcout << "GMWM estimator: " << GMWM << std::endl;
  
  return GMWM;                          
}

//' @title GMWM Step
//' @details
//' Allows the user to supply their "guess" instead of a random draw.
//' @example
//' x=rnorm(100)
// [[Rcpp::export]]
arma::mat GMWM_adv(const arma::vec& theta, const std::vector<std::string>& desc,
                   const arma::mat& V, const arma::vec& wv_empir,
                   const arma::vec& tau, unsigned int N){
                                 
       // Number of parameters
  unsigned int num_param = theta.n_elem;
  
  // Initialisation of results structures
  arma::mat GMWM(1,num_param);
    
  // Starting values
  arma::vec starting_theta = set_starting_values(theta, desc);
  
  starting_theta = Rcpp_OptimStart(starting_theta, desc, tau, wv_empir, N);
      
  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
      
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, omega, tau, wv_empir, N);

  // Save results
  GMWM.row(0) = set_result_values(estim_GMWM, desc);
  
  return GMWM;                          
}



//' @title Simulate GMWM
//' 
//' @example
//' x=rnorm(100)
//' wavelet_variance_arma(x, "haar", "diag")
// [[Rcpp::export]]
arma::field<arma::mat> simGMWM(const arma::vec& theta, const arma::mat& omega,
                               const std::vector<std::string>& desc, const arma::vec& wv_empir,
                               const arma::vec& tau, unsigned int N, unsigned int B = 500, bool var_or_mu = false){
  
  // Number of parameters
  unsigned int num_param = theta.n_elem;
  
  // Initialisation of results structures
  arma::mat GMWM(B,num_param);
  arma::mat GMWM_plus(B,num_param);
  
  // Starting values
	arma::vec starting_theta = set_result_values(theta, desc);

  // Start bootstrap
  for(unsigned int b=0; b<B; b++){  	
  	// Simulate  	
  	arma::vec x = gen_model(N, theta, desc);
    
  	// ------------------------------------
  	// Compute standard GMWM
  	// ------------------------------------
  	// Compute WV
  	arma::field<arma::mat> wv_x = wavelet_variance_arma(x, "haar", "diag");
  	
  	// Omega matrix
  	arma::mat omega = arma::inv(diagmat(wv_x(4)));
  	
  	// Empirical WV
  	arma::vec wv_empir = wv_x(0);
      	
  	// Find GMWM estimator
  	arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, omega, tau, wv_empir, N);
  	
  	// Save results
  	GMWM.row(b) = set_result_values(estim_GMWM, desc);
  	
  	// ------------------------------------
  	// Compute augmented GMWM
  	// ------------------------------------
  	// Compute Omega
  	arma::mat V = gmwm_bootstrapper(GMWM.row(b), desc, max(tau), N, B, var_or_mu);
  	arma::mat omega_v = arma::inv(diagmat(V));
  	
  	// Empirical WV + variance
    arma::vec temp(1);
    if(var_or_mu){
      temp(0) = arma::var(x);
    }
    else{
      temp(0) = arma::mean(x);
    }
    
  	arma::vec wv_empir_v = join_cols(temp,wv_empir);
  	
  	// Find GMWM+ estimator
  	arma::vec estim_GMWM_plus = Rcpp_Optim(theta, desc, omega_v, tau, wv_empir_v, N);
  	
  	// Save results
  	GMWM_plus.row(b) = set_result_values(estim_GMWM_plus, desc);
  }
  
  arma::field<arma::mat> out(2);
  out(0) = GMWM;
  out(1) = GMWM_plus;
  
  return out;
} 