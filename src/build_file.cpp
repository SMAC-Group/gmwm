#include <RcppArmadillo.h>
#include <string>
#include <map>
// #include <omp.h>

using namespace arma;
using namespace Rcpp;

// // [[Rcpp::plugins(openmp)]]

// Support functions

#define my_isok(x) (!ISNA(x) & !ISNAN(x))

inline double square(double x){
  return x*x;
}

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
//' @param x A \code{vec} containing probabilities (e.g. 0 <= x <= 1)
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




//' @title Reverse Subset Column
//' @description Subsets the column by going from high indices to low (the reverse of the supported practice)
//' @usage rev_col_subset(x, start, end)
//' @param x A \code{matrix} of dimensions M x N
//' @param start A \code{unsigned int} that indicates the starting column.
//' @param end A \code{unsigned int} that indicates the ending column.
//' @return x A \code{matrix} with matrix rows displayed in reverse order
//' @details Consider a vector x=[[1,2],[3,4]].
//' By setting \code{start=1} and \code{end=0}, the function would output x=[[2,1],[4,1]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix cols start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow=2,byrow=T)
//' rev_col_subset(x, 1, 0)
// [[Rcpp::export]]
arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(x.n_rows, start-end+1);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.col(i) = x.col(start-i);
  }
  return A;
}

// [[Rcpp::export]]
arma::vec diff_cpp(const arma::vec& x, unsigned int lag = 1){
  unsigned int n=x.n_elem;
  return (x.rows(lag,n-1) - x.rows(0,n-lag-1));
}


//' @title Reverse Subset Row
//' @description Subsets the row by going from high indices to low (the reverse of the supported practice)
//' @usage rev_row_subset(x, start, end)
//' @param x A \code{matrix} of dimensions M x N
//' @param start A \code{unsigned int} that indicates the starting row.
//' @param end A \code{unsigned int} that indicates the ending row.
//' @return x A \code{matrix} with matrix rows displayed in reversed order
//' @details Consider a vector x=[[1,2],[3,4]], the function would output x=[[3,4],[1,2]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix rows start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow=2,byrow=T)
//' rev_row_subset(x, 1, 0)
// [[Rcpp::export]]
arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(start-end+1, x.n_cols);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.row(i) = x.row(start-i);
  }
  return A;
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


//' @title Converting an ARMA Process to an Infinite MA Process
//' @description Takes an ARMA function and converts it to an infinite MA process.
//' @usage ARMAtoMA_cpp(ar, ma, lag_max)
//' @param ar A \code{column vector} of length p
//' @param ma A \code{column vector} of length q
//' @param lag_max A \code{int} of the largest MA(Inf) coefficient required.
//' @return A \code{column vector} containing coefficients
//' @details This function is a port of the base stats package's ARMAtoMA. There is no significant speed difference between the two.
//' @author R Core Team and JJB
//' @examples
//' # ARMA(2,1)
//' ARMAtoMA_cpp(c(1.0, -0.25), 1.0, 10)
//' # ARMA(0,1)
//' ARMAtoMA_cpp(numeric(0), 1.0, 10)
// [[Rcpp::export]]
arma::vec ARMAtoMA_cpp(arma::vec ar, arma::vec ma, int lag_max)
{
  int p = ar.n_elem;
  int q = ma.n_elem;
  int m = lag_max;
  
  double tmp;
  
  arma::vec psi(m);
  
  if(m <= 0 || m == NA_INTEGER){
    Rcpp::stop("invalid value of lag.max");
  }
  
  for(int i = 0; i < m; i++) {
    tmp = (i < q) ? ma(i) : 0.0;
    for(int j = 0; j < std::min(i+1, p); j++){
      tmp += ar(j) * ((i-j-1 >= 0) ? psi(i-j-1) : 1.0);
    }
    psi(i) = tmp;
  }
  return psi;
}

//' @title Time Series Convolution Filters
//' @description Applies a convolution filter to a univariate time series.
//' @usage cfilter(x, filter, sides, circular)
//' @param x A \code{column vector} of length T
//' @param filter A \code{column vector} of length f
//' @param sides An \code{int} that takes either 1:for using past values only or 2: filter coefficients are centered around lag 0.
//' @param circular A \code{bool} that indicates if the filter should be wrapped around the ends of the time series.
//' @return A \code{column vec} that contains the results of the filtering process.
//' @details This is a port of the cfilter function harnessed by the filter function in stats. 
//' It is about 5-7 times faster than R's base function. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
//' @author R Core Team and JJB
//' @examples
//' x = 1:100
//' # 
//' cfilter(x, rep(1, 3))
//' # Using R's function
//' filter(x, rep(1, 3))
//' #
//' cfilter(x, rep(1, 3), sides = 1)
//' # Using R's function
//' filter(x, rep(1, 3), sides = 1)
//' #
//' cfilter(x, rep(1, 3), sides = 1, circular = TRUE)
//' # Using R's function
//' filter(x, rep(1, 3), sides = 1, circular = TRUE)
// [[Rcpp::export]]
arma::vec cfilter(arma::vec x, arma::vec filter, int sides = 2, bool circular = false)
{
  
  int nx = x.n_elem;
  int nf = filter.n_elem;
  int nshift;
  
  if(sides == NA_INTEGER || circular == NA_LOGICAL)  Rcpp::stop("invalid input");
  
  double z, tmp;
  
  if(sides == 2){
    nshift = nf /2;
  }
  else{
    nshift = 0;
  }
    
  arma::vec out = arma::zeros<arma::vec>(nx);
  
  if(!circular) {
    for(int i = 0; i < nx; i++) {
      z = 0;
      if(i + nshift - (nf - 1) < 0 || i + nshift >= nx) {
        out(i) = NA_REAL;
        continue;
      }
      for(int j = std::max(0, nshift + i - nx); j < std::min(nf, i + nshift + 1) ; j++) {
        tmp = x(i + nshift - j);
        z += filter(j) * tmp;
      }
      out(i) = z;
    }
  } else { /* circular */
  for(int i = 0; i < nx; i++)
  {
    z = 0;
    for(int j = 0; j < nf; j++) {
      int ii = i + nshift - j;
      if(ii < 0) ii += nx;
      if(ii >= nx) ii -= nx;
      tmp = x(ii);
      z += filter(j) * tmp;
    }
    out(i) = z;
  }
  }
  return out;
}


//' @title Time Series Recursive Filters
//' @description Applies a recursive filter to a univariate time series.
//' @usage rfilter(x, filter, init)
//' @param x A \code{column vector} of length T
//' @param filter A \code{column vector} of length f
//' @param init A \code{column vector} of length f that contains the initial values of the time series in reverse.
//' @return x A \code{column vector} with its contents reversed.
//' @details Note: The length of 'init' must be equal to the length of 'filter'.
//' This is a port of the rfilter function harnessed by the filter function in stats. 
//' It is about 6-7 times faster than R's base function. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
//' @author R Core Team and JJB
//' @examples
//' x = 1:100
//' # 
//' rfilter(x, rep(1, 3), rep(1, 3))
//' # Using R's function
//' filter(x, rep(1, 3), method="recursive", init=rep(1, 3))
// [[Rcpp::export]]
arma::vec rfilter(arma::vec x, arma::vec filter, arma::vec init)
{
 
  int nx = x.n_elem, nf = filter.n_elem;
    
  double sum;
  arma::vec r = arma::join_cols(reverse_vec(init), arma::zeros<arma::vec>(nx) ); 
  // see filter source
  // .Call(C_rfilter, x, filter, c(rev(init[, 1L]), double(n)))[-ind]
  // r is then c(rev(init[, 1L]), double(n))
  arma::vec rx = x;
  arma::vec rf = filter;
  
  for(int i = 0; i < nx; i++) {
    sum = rx(i);
    for (int j = 0; j < nf; j++) {
      if(nf + i - j - 1 >= 0){
        sum += r(nf + i - j - 1) * rf(j);
      }else{
        r[nf + i] = NA_REAL; goto bad3; 
      }
    }
    r(nf + i) = sum;
    bad3:
  continue;
  }

  return r.rows(nf,r.n_elem-1); // returns truncated vec (discards filter)
}

//' @title Expand Grid for Same Dimensional Case
//' @description Creates the different pairings possible with two different variables.
//' @usage expand_grid_red(nx)
//' @param nx An \code{integer} of length f that contains the initial values of the time series in reverse.
//' @return x A \code{matrix} listing values from 1...nx in one column and 1...1, 2...2,....,n...n, in the other
//' @author JJB
//' @details This function is hidden and is not accessible from R.
//' @name expand_grid_red
//' @docType methods
//' @rdname expand_grid_red-methods
//' @keywords internal
arma::mat expand_grid_red(int nx){
  
  arma::mat g(nx*nx,2);
  int j = 1;
  
  for(int i = 1; i <= nx*nx; i++){
    int mod = i%nx;
    if(mod == 0){
      mod = nx;
    }
    g(i-1,0) = mod;
    g(i-1,1) = j;
    if(i%nx == 0){
      j++;
    }
  }
  
  return g;
}


//' @title Compute Theoretical ACF for an ARMA Process
//' @description Compute the theoretical autocorrelation function for an ARMA process.
//' @usage ARMAacf_cpp(ar,ma,lag_max)
//' @param ar A \code{vector} of length p containing AR coefficients
//' @param ma A \code{vector} of length q containing MA coefficients
//' @param lag_max A \code{unsigned integer} indicating the maximum lag necessary
//' @return x A \code{matrix} listing values from 1...nx in one column and 1...1, 2...2,....,n...n, in the other
//' @details This is an implementaiton of the ARMAacf function in R. It is approximately 40x times faster. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
//' @author R Core Team and JJB
//' @examples
//' # ARMA(2,1)
//' ARMAacf_cpp(c(1.0, -0.25), 1.0, lag_max = 10)
//' # ARMA(0,1)
//' ARMAacf_cpp(numeric(0), .35, lag_max = 10)
// [[Rcpp::export]]
arma::vec ARMAacf_cpp(arma::vec ar, arma::vec ma, unsigned int lag_max) 
{
  unsigned int p = ar.n_elem;
  unsigned int q = ma.n_elem;
  
  arma::vec Acf;
  if (!p && !q){
    Rcpp::stop("empty model supplied");
  }
  unsigned int r = std::max(p, q + 1);
  if (p > 0) {
    if (r > 1) {
      if (r > p) {
        unsigned int temp_size = ar.n_elem;
        ar.resize(temp_size+r-p);
        p = r;
      }
      
      arma::mat A = arma::zeros<arma::mat>(p + 1, 2 * p + 1);
            
      // A[ind] <- c(1, -ar)
      arma::rowvec temp = arma::rowvec(p + 1);
      temp(0) = 1; 
      temp.cols(1,p) = -1*arma::conv_to<arma::rowvec>::from(ar);
            
      for(unsigned int i = 0; i <= p; i++){
        //row,col, row, col
        A.submat(i,i,i,i+p) = temp;
      }
            
      // A[, (2L * p + 1L):(p + 2L)]
      arma::mat tempmat = rev_col_subset(A,2 * p,p + 1); // start to end
        
      // A[, 1L:p] <- A[, 1L:p] + A[, (2L * p + 1L):(p + 2L)]
      A.cols(0,p-1) = A.cols(0,p-1) + tempmat;
            
      // rhs <- c(1, rep(0, p))
      arma::vec rhs = arma::zeros<arma::vec>(p+1);
      rhs(0) = 1;
            
      if (q > 0) {
        arma::vec psi = arma::vec(q+1);
        psi(0) = 1;
        psi.rows(1,q) = ARMAtoMA_cpp(ar, ma, q);
                
        arma::vec theta = arma::zeros<arma::vec>(2*q+2);
        theta(0) = 1;
        theta.rows(1,q) = ma;

        // in 1 + 0:q
        for (unsigned int k = 0; k <= q; k++){
          rhs(k) = sum(psi % theta.rows(k,k+q));
        }
      }
      Acf = solve(rev_row_subset(rev_col_subset(A,p,0),p,0), rhs);
      Acf = Acf.rows(1,Acf.n_rows-1)/Acf(0);
    }
    else{
      Acf = ar;
    }
    
    if (lag_max > p) {
      arma::vec xx = arma::zeros<arma::vec>(lag_max - p);
      Acf = arma::join_cols(Acf, rfilter(xx, ar, reverse_vec(Acf) ) );
    }
    
    //  Acf = c(1, Acf[1L:lag.max])
    arma::vec temp = arma::vec(lag_max+1);
    temp(0) = 1;
    temp.rows(1,lag_max) = Acf.rows(0,lag_max-1);
    
    Acf = temp;
  }
  else if (q > 0) {
    
    //  x = c(1, ma)
    arma::vec x(q+1);
    x(0) = 1;
    x.rows(1,q) = ma;

    Acf = cfilter(arma::join_cols(x, arma::zeros<arma::vec>(q)), reverse_vec(x), 1);
    
    // [-(1L:q)]
    Acf = Acf.rows(q,Acf.n_elem-1);
    
    if (lag_max > q){
      Acf = arma::join_cols(Acf, arma::zeros<arma::vec>(lag_max - q));
    }
    Acf = Acf/Acf(0);
  }
  
  return Acf;
}

//' @title ARMA process to WV
//' @description This function computes the (haar) WV of an ARMA process
//' @param ar A \code{vec} containing the coefficients of the AR process
//' @param ma A \code{vec} containing the coefficients of the MA process
//' @param tau A \code{vec} containing the scales e.g. 2^tau
//' @param sigma A \code{double} containing the residual variance
//' @return A \code{vec} containing the wavelet variance of the ARMA process.
//' @examples
//' arma_to_wv(c(.23,.43), c(.34,.41,.59), 2^(1:9), 3)
//' @seealso \code{\link{ARMAtoMA_cpp}},\code{\link{ARMAacf_cpp}}
// [[Rcpp::export]]
arma::vec arma_to_wv(arma::vec ar, arma::vec ma, arma::vec tau, double sigma) {
  
  arma::vec n = arma::sort(tau/2);
  unsigned int ntau = tau.n_elem;
  double sig2 = (arma::sum(arma::square(ARMAtoMA_cpp(ar,ma,1000)))+1)*sigma;
  
  arma::vec wvar(ntau);
  
  // initial starting term
  arma::vec term4 = ARMAacf_cpp(ar, ma, n(0));
  wvar(0)=( ( ( n(0)*(1.0-term4(term4.n_elem-1))) / square(n(0)))*sig2)/2.0;

  for(unsigned int j = 1; j < ntau; j++){
    arma::vec boh(n(j) - 1);
    for (int i=1; i<= n(j) - 1; i++){
      arma::vec term1=ARMAacf_cpp(ar, ma, (n(j)-i));
      arma::vec term2=ARMAacf_cpp(ar, ma, i);
      arma::vec term3=ARMAacf_cpp(ar, ma, (2*n(j)-i));
      // Account for starting loop at 1 instead of 0.
      boh(i-1)=i*((2.0*term1(term1.n_elem-1))-term2(term2.n_elem-1)-term3(term3.n_elem-1));
    }
    arma::vec term4=ARMAacf_cpp(ar, ma, n(j));
    wvar(j)=((( (n(j)*(1.0-term4(term4.n_elem-1)) ) + arma::sum(boh) ) /square(n(j)))*sig2)/2.0;
  }
  
  return wvar;
}


//' @title Compute Tau-Overlap Allan Variance
//' @description Computation of Tau-Overlap Allan Variance
//' @usage avar_to_cpp(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{matrix} that contains:
//' \itemize{
//'  \item{Col 1}{The size of the cluster}
//'  \item{Col 2}{The Allan variance}
//'  \item{Col 3}{The error associated with the variance estimation.}
//' }
//' @details
//' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
//' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
//' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
//' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist. 
//' The tau-overlap estimator is given by:
//' 
//' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
//' 
//' @author JJB
//' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avar_to_cpp(combined.ts)
// [[Rcpp::export]]
arma::mat avar_to_cpp(arma::vec x) {
  
   // Length of vector
   unsigned int T = x.n_elem;
   // Create the number of halves possible and use it to find the number of clusters
   unsigned int J = floor(log10(T)/log10(2))-1;
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 1; i <= J; i++){
     // Tau
     unsigned int tau = pow(2,i);

     // Y.Bar
     unsigned int N = floor(T/tau);
     arma::vec yBar = arma::zeros<arma::vec>(N);
     for(unsigned int j = 0; j < N;j++){
       yBar(j) = sum( x.rows(tau*j, tau*j+tau - 1) )/tau;
     }

     // Clusters
     unsigned int M = floor(T/(2*tau) );
  	 double summed = 0;
		 for(unsigned int k = 0; k < M; k++){
			 summed +=  pow(yBar(2*k+1) - yBar(2*k),2);
		 }
     
     // Cluster size
     av(i-1,0) = tau; 
     // Compute the Allan Variance estimate
     av(i-1,1) = summed/(2*M); 
     // Compute Error
     av(i-1,2) = 1/sqrt(2*( (double(T)/tau) - 1) );
  }
  
  // Return as list
  return av;
}

//' @title Compute Maximal-Overlap Allan Variance using Means
//' @description Computation of Maximal-Overlap Allan Variance
//' @usage avar_mo_cpp(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"allan"}{The Allan variance}
//'  \item{"errors"}{The error associated with the variance estimation.}
//' }
//' @details
//' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
//' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
//' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
//' Then, \eqn{M = N - 2n} samples exist. 
//' The Maximal-overlap estimator is given by:
//' \eqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }
//' 
//' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
//' @author JJB
//' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avar_mo_cpp(combined.ts)
// [[Rcpp::export]]
arma::mat avar_mo_cpp(arma::vec x) {
  
   // Length of vector
   unsigned int T = x.n_elem;
   
   // Create the number of halves possible and use it to find the number of clusters
   unsigned int J = floor(log10(T)/log10(2))-1;
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 1; i <= J; i++){
     // Tau
     unsigned int tau = pow(2,i);

     // Y.Bar
     arma::vec yBar = arma::zeros<arma::vec>(T);
     for(unsigned int j = 0; j <= T-tau; j++){
       yBar(j) = sum( x.rows(j, tau+j -1) ) / tau;
     }
     
     // Clusters
     unsigned int M = T-2*tau;
		 double summed = 0;
		 for(unsigned int k = 0; k < M; k++){
			 summed +=  pow(yBar(k) - yBar(k+tau),2);
		 }
     
     // Cluster size
     av(i-1,0) = tau; 
     // Compute the Allan Variance estimate
     av(i-1,1) = summed/(2*(T - 2*tau + 1)); 
     // Compute Error
     av(i-1,2) = 1/sqrt(2*( (double(T)/tau) - 1) );
  }
  
  return av;
}

//' @title Generate a white noise process
//' @description Generates a white noise process with variance parameter sigma.
//' @param N An \code{integer} for signal length.
//' @param sigma_WN A \code{double} that contains process standard deviation.
//' @return wn A \code{vec} containing the white noise.
//' @examples
//' gen_wn(10, 1.5)
// [[Rcpp::export]]
arma::vec gen_wn(const unsigned int N, const double sigma_WN)
{
	arma::vec wn(N);
  
  for(unsigned int i = 0; i < N; i++){
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
//' @param sigma_rw A \code{double} for the deviation of the random walk.
//' @return grw A \code{vec} containing the random walk without drift.
//' @examples
//' gen_rw(10, 8.2)
// [[Rcpp::export]]
arma::vec gen_rw(const unsigned int N, const double sigma_rw = 1)
{
  arma::vec grw(N);
  for(unsigned int i = 0; i < N; i++){
      grw(i) = R::rnorm(0.0, sigma_rw);
  }
  return cumsum(grw);
}

//' @title Quantisation Noise to WV
//' @description This function compute the WV (haar) of a Quantisation Noise (QN) process
//' @param q2 A \code{double} corresponding to variance of drift
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the QN.
//' @examples
//' x.sim = 1:1000
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = qn_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec qn_to_wv(double q2, const arma::vec& Tau){
  return 3*q2/(2*arma::square(Tau));
}

//' @title White Noise to WV
//' @description This function compute the WV (haar) of a White Noise process
//' @param sig2 A \code{double} corresponding to variance of WN
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the white noise.
//' @examples
//' x.sim = cumsum(rnorm(100000))
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = wn_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec wn_to_wv(double sig2, arma::vec Tau){
	return sig2/Tau;
}


//' @title Random Walk to WV
//' @description This function compute the WV (haar) of a Random Walk process
//' @param sig2 A \code{double} corresponding to variance of RW
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the random walk.
//' @examples
//' x.sim = cumsum(rnorm(100000))
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = rw_to_wv(1,tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec rw_to_wv(double sig2, const arma::vec& Tau){
	return sig2*((2*arma::square(Tau) + 2)/(24*Tau));
}


//' @title Drift to WV
//' @description This function compute the WV (haar) of a Drift process
//' @param omega A \code{double} corresponding to variance of drift
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the drift.
//' @examples
//' x.sim = 1:1000
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = dr_to_wv(1, tau)
//' plot(tau, wv.theo, col = "red")
// [[Rcpp::export]]
arma::vec dr_to_wv(double omega,const arma::vec& Tau){
	return square(omega)*arma::square(Tau)/16;
}

//' @title AR1 process to WV
//' @description This function compute the WV (haar) of an AR(1) process
//' @param phi A \code{double} that is the phi term of the AR(1) process
//' @param sig2 A \code{double} corresponding to variance of AR(1) process
//' @param Tau A \code{vec} containing the scales e.g. 2^tau
//' @return A \code{vec} containing the wavelet variance of the AR(1) process.
//' @examples
//' x.sim = gen_ar1( N = 10000, phi = 0.9, sig2 = 4 )
//' ntau = floor(log(length(x.sim),2))
//' tau = 2^(1:ntau)
//' wv.theo = ar1_to_wv(phi = 0.9, sig2 = 16, tau)
//' plot(tau, wv.theo, col = "red")
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
//' @examples
//' e_drift(1,200)
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
//' @examples
//' m2_drift(1, 200)
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
//' @examples
//' var_drift(1, 200)
// [[Rcpp::export]]
arma::vec var_drift(double omega, int n_ts){
  // Compute m1
	arma::vec m1 = e_drift(omega, n_ts);
	
	// Compute m2
	arma::vec m2 = m2_drift(omega, n_ts);
	
	// Compute var
  return (m2 - m1*m1)*double(n_ts)/double(n_ts-1.0);
}


//' @title Quadrature Mirror Filter
//' @description Calculate the series quadrature mirror filter (QMF). Requires a series of an even length.
//' @usage qmf(g, inverse)
//' @param g A \code{vector} that contains the filter constants.
//' @param inverse A \code{bool} that indicates whether the inverse quadrature mirror filter is computed. 
//' By default, the inverse quadrature mirror is computed.
//' @return A \code{vector} that contains either the forward QMF (evalute in order) or the inverse QMF (reverse order). 
//' @author JJB
//' @examples
//' # Haar values
//' g = rep(1/sqrt(2),2)
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
//' This template can be used to increase the amount of filters available for selection.
//' @author JJB
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
//' @description Calculation of the coefficients for the discrete wavelet transformation. 
//' @usage dwt_cpp(x, filter_name, nlevels, boundary)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param filter_name A \code{string} indicating the filter.
//' @param nlevels An \code{integer} indicating the level of the decomposition.
//' @param boundary A \code{string} indicating the type of boundary method to use. Either \code{boundary="periodic"} or \code{"reflection"}.
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm
//' @author JJB
//' @examples
//' set.seed(999)
//' x = rnorm(100)
//' dwt_cpp(x, "haar", 4, boundary="periodic")
// [[Rcpp::export]]
arma::field<arma::vec> dwt_cpp(arma::vec x, String filter_name = "haar", 
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
arma::field<arma::vec> modwt_cpp(arma::vec x, String filter_name = "haar", 
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
//' The modulus is defined as: \eqn{r = Mod\left(z\right) = \sqrt{\left(x^2 + y^2\right)}}{r = Mod(z) = sqrt(x^2 + y^2)}
//' This function will return: \eqn{r^2 = Mod\left(z\right)^2 = x^2 + y^2}
//' @examples
//' Mod_squared_cpp(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_squared_cpp(arma::cx_vec x){
   return arma::square(arma::real(x)) + arma::square(arma::imag(x));
}

//' @title Absolute Value or Modulus of a Complex Number.
//' @description Computes the value of the Modulus.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the modulus for each element.
//' @details Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
//' The modulus is defined as: \eqn{r = Mod(z) = \sqrt{(x^2 + y^2)}}
//' @examples
//' Mod_cpp(c(1+.5i, 2+1i, 5+9i))
// [[Rcpp::export]]
arma::vec Mod_cpp(arma::cx_vec x){
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
//' @examples
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
//' The vectors are truncated by removing the first n wavelet coefficients. 
//' These vectors are then stored into the field that is returned.
//' Note: As a result, there are no NA's introduced and hence the na.omit is not needed.
//' @examples
//' x=rnorm(100)
//' brick_wall(modwt_cpp(x))
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

// [[Rcpp::export]]
arma::mat field_to_matrix(arma::field<arma::vec> x, unsigned int row){
  unsigned int nx = x.n_elem;
  arma::mat A(row,nx);
  for(unsigned int i =0; i<nx; i++){
    A.col(i) = x(i);
  }
  return A; 
}


// Implemented for Robust

// @title Objective Function for Tuning Constant
// @description Objective function that finds tuning constant
// @usage objFun_find_biwc(crob, eff)
// @param crob A \code{double} that indicates the value to optimize. Normally between [0,5].
// @param eff A \code{double} that represents the efficiency.
// @return A \code{double} that contains the biwc.
// @details This function is not available from within R.
// @name objFun_find_biwc
// @docType methods
// @rdname objFun_find_biwc-methods
double objFun_find_biwc(double crob, double eff){
  
  // q, mean, sd, lower.tail, log.p
  double norm_prob = R::pnorm(crob,0.0,1.0,1,0);
  double norm_density = R::dnorm(crob,0.0,1.0,0);
  
  // avoid using pow as much as possible
  double crob_sq = square(crob);
  double crob_quad = square(crob_sq);
  double crob_eight = square(crob_quad);
  
  double mu20c=1309458150*norm_prob-2*crob*(654729075+218243025*crob_sq+43648605*crob_quad+6235515*(crob_sq*crob_quad)+692835*crob_eight+62985*(crob_eight*crob_sq)+4845*(crob_quad*crob_eight)+323*(crob_sq*crob_quad*crob_eight)+19*(crob_eight*crob_eight)+(crob_eight*crob_eight*crob_sq))*norm_density-654729075;
  double mu18c=68918850*norm_prob-2*crob*(34459425+11486475*crob_sq+2297295*crob_quad+328185*(crob_sq*crob_quad)+36465*crob_eight+3315*(crob_eight*crob_sq)+255*(crob_quad*crob_eight)+17*(crob_sq*crob_quad*crob_eight)+(crob_eight*crob_eight))*norm_density-34459425;
  double mu16c=4054050*norm_prob-2*crob*(2027025+675675*crob_sq+135135*crob_quad+19305*(crob_sq*crob_quad)+2145*crob_eight+195*(crob_eight*crob_sq)+15*(crob_quad*crob_eight)+(crob_sq*crob_quad*crob_eight))*norm_density-2027025;
  double mu14c=270270*norm_prob-2*crob*(135135+45045*crob_sq+9009*crob_quad+1287*(crob_sq*crob_quad)+143*crob_eight+13*(crob_eight*crob_sq)+(crob_quad*crob_eight))*norm_density-135135;
  double mu12c=20790*norm_prob-2*crob*(10395+3465*crob_sq+693*crob_quad+99*(crob_sq*crob_quad)+11*crob_eight+(crob_eight*crob_sq))*norm_density-10395;
  double mu10c=1890*norm_prob-2*crob*(945+315*crob_sq+63*crob_quad+9*(crob_sq*crob_quad)+crob_eight)*norm_density-945;
  double mu8c=210*norm_prob-2*crob*(105+35*crob_sq+7*crob_quad+(crob_sq*crob_quad))*norm_density-105;
  double mu6c=30*norm_prob-2*crob*(15+5*crob_sq+crob_quad)*norm_density-15;
  double mu4c=6*norm_prob-2*crob*(3+crob_sq)*norm_density-3;
  double mu2c=2*norm_prob-2*crob*norm_density-1;
  double ac=(1/crob_eight)*mu10c-(4/(crob_sq*crob_quad))*mu8c+(6/crob_quad)*mu6c-(4/crob_sq)*mu4c+mu2c;
  double Q=(28/crob_quad)*mu8c-(8/crob_sq)*mu6c+(1/(crob_eight*crob_eight))*mu20c+(8/(crob_sq*crob_quad*crob_eight))*mu18c+(28/(crob_quad*crob_eight))*mu16c-(56/(crob_eight*crob_sq))*mu14c+(70/crob_eight)*mu12c-(56/(crob_sq*crob_quad))*mu10c+mu4c-ac*ac;
  double M=(1/crob_eight)*mu12c-(4/(crob_sq*crob_quad))*mu10c+(6/crob_quad)*mu8c-(4/crob_sq)*mu6c+mu4c-ac;
  return square((0.5*M*M)/Q - eff);
}

// @title Obtain Tuning Constant crob.bw
// @description Objective function that finds tuning constant
// @usage find_biwc(eff)
// @param eff A \code{double} that represents the desired level of efficiency as compared to the classic MLE.
// @return A \code{double} that contains the optimized biwc.
// @examples
// find_biwc(0.6)
double find_biwc(double eff){
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optimize = stats["optimize"];    

  Rcpp::List crob_t = optimize(_["f"]  = Rcpp::InternalFunction(&objFun_find_biwc),
                          _["lower"] = 0,
                          _["upper"] = 10,
                          _["eff"] = eff);
                          
  double out = as<double>(crob_t[0]);

  return out;
}


// Objective funcion for robust WV
double objFun_sig_rob_bw(double sig2_bw, arma::vec x, double a_of_c, double crob_bw){
  arma::vec r = x/sqrt(sig2_bw);
  arma::vec rsq = arma::square(r);
  arma::uvec ivec = (abs(r) > crob_bw);
  arma::vec w = ((1 - ivec) % arma::square(1 - rsq/square(crob_bw)) );
  return square(arma::mean(arma::square(r)%arma::square(w)) - a_of_c);
}

// Robust estimator. Inputs are wavelet coefficients (y) and desired level of efficiency. 
double sig_rob_bw(arma::vec y, double eff=0.6){
  double crob_bw = find_biwc(eff);

  arma::vec x = y/arma::stddev(y);  
  
  // q, mean, sd, lower.tail, log.p
  double norm_prob = R::pnorm(crob_bw,0.0,1.0,1,0);
  double norm_density = R::dnorm(crob_bw,0.0,1.0,0);
  
  // avoid using pow as much as possible
  double crob_sq = square(crob_bw);
  double crob_quad = square(crob_sq);
  double crob_eight = square(crob_quad);
  
  double a_of_c = (1/crob_eight)*(1890*norm_prob-2*crob_bw*(945+315*crob_sq+63*crob_quad+9*(crob_sq*crob_quad)+crob_eight)*norm_density-945)
                  -(4/(crob_sq*crob_quad))*(210*norm_prob-2*crob_bw*(105+35*crob_sq+7*crob_quad+(crob_sq*crob_quad))*norm_density-105)
                  +(6/crob_quad)*(30*norm_prob-2*crob_bw*(15+5*crob_sq+crob_quad)*norm_density-15)
                  -(4/crob_sq)*(6*norm_prob-2*crob_bw*(3+crob_sq)*norm_density-3)
                  +2*norm_prob-2*crob_bw*norm_density-1;

  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optimize = stats["optimize"];    

  Rcpp::List opt_val = optimize(_["f"]  = Rcpp::InternalFunction(&objFun_sig_rob_bw),
                          _["lower"] = 0,
                          _["upper"] = 5,
                          _["x"] = x,
                          _["a_of_c"] = a_of_c,
                          _["crob_bw"] = crob_bw);
          
  double sig2_hat_rob_bw = as<double>(opt_val[0])*var(y);
  return sig2_hat_rob_bw;
}

//' @title Generate eta3 confidence interval
//' @description Computes the eta3 CI
//' @param y A \code{vec} that computes the brickwalled modwt dot product of each wavelet coefficient divided by their length.
//' @param dims A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @examples
//' x=rnorm(100)
//' out = brick_wall(modwt_cpp(x), haar_filter())
//' sample = out[[1]]
//' y = sample*sample/length(sample)
//' ci_eta3(y, length(sample), 0.025)
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

// [[Rcpp::export]]
arma::mat ci_eta3_robust(arma::vec y, arma::vec dims, double p, double eff) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);
    eff = sqrt(eff);
    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2,i+1),1.0);
      out(i,1) = eff * eta3 * y(i)/(R::qchisq(1-p, eta3, 1, 0)); // Lower CI
      out(i,2) = eta3 * y(i)/(eff*R::qchisq(p, eta3, 1, 0)); // Upper CI
    }

    out.col(0) = y;

    return out;
}

//' @title Generate a Confidence intervval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a chi-squared confidence interval
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param y A \code{vec} that contains the wave variance.
//' @param type A \code{String} indicating the confidence interval being calculated.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level.
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the structure:
//' \itemize{
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' }
//' @details 
//' This function can be expanded to allow for other confidence interval calculations.
//' @examples
//' x=rnorm(100)
//' wave_variance(brick_wall(modwt_cpp(x), haar_filter()))
// [[Rcpp::export]]
arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& y, String type = "eta3", double p = 0.025, bool robust = false, double eff = 0.6){
    
  unsigned int nb_level = y.n_elem;
  arma::vec dims(nb_level);
  
  for(unsigned int i = 0; i < nb_level; i++){
    dims(i) = signal_modwt_bw(i).n_elem;
  }

  arma::mat out(y.n_elem , 3);

  
  if(type == "eta3"){
      if(robust){
        out = ci_eta3_robust(y, dims, p, eff);
      }
      else{
        out = ci_eta3(y, dims,p);  
      }
  }
  else{
      stop("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}


//' @title Generate a Wave Variance for a Univariate Time Series
//' @description Computes an estimate of the wave variance
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param robust A \code{boolean} to determine the type of wave estimation.
//' @param eff A \code{double} that indicates the efficiency.
//' @return A \code{vec} that contains the wave variance.
//' @examples
//' x=rnorm(100)
//' wave_variance(brick_wall(modwt_cpp(x), haar_filter()))
// [[Rcpp::export]]
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw, bool robust = false, double eff = 0.6){
  
  unsigned int nb_level = signal_modwt_bw.n_elem;
  arma::vec y(nb_level);
  
  if(robust){
    // Robust wavelet variance estimation
    for(unsigned int i=0; i < nb_level; i++){
      arma::vec wav_coef = sort(signal_modwt_bw(i));
      y(i) = sig_rob_bw(wav_coef, eff);
    }
  }else{
    // Classical wavelet variance estimation
    for(unsigned int i=0; i < nb_level;i++){
      arma::vec temp = signal_modwt_bw(i);
      y(i) = dot(temp,temp)/temp.n_elem;
    }
  }
  
  return y;
}


//' @title Computes the (MODWT) wavelet variance
//' @description Calculates the (MODWT) wavelet variance
//' @param signal_modwt A \code{field<vec>} that contains the modwt decomposition.
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @param robust A \code{boolean} that triggers the use of the robust estimate.
//' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
//' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @param ci_type A \code{String} indicating the confidence interval being calculated. Valid value: "eta3"
//' @param strWavelet A \code{String} indicating the type of wave filter to be applied. Must be "haar"
//' @return A \code{mat} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance},
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @details 
//' This function powers the wvar object. It is also extendable...
//' @examples
//' x=rnorm(100)
//' decomp = modwt(x)
//' wvar_cpp(decomp$data, decomp$nlevels, robust=false, eff=0.6, p = 0.025, strWavelet = "haar")
// [[Rcpp::export]]
arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt, unsigned int nb_level, bool robust=false, double eff=0.6, double p = 0.025, std::string ci_type="eta3", std::string strWavelet="haar") {
  
  // MODWT transform
  arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, select_filter(strWavelet));
  
  // Wavelet Variance
  arma::vec y = wave_variance(signal_modwt_bw, robust, eff);
  
  // Confidence Interval
  return ci_wave_variance(signal_modwt_bw, y, ci_type, p, robust, eff);
}

//' @title Computes the MODWT scales
//' @description Calculates the MODWT scales
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @return A \code{vec} that contains 2^1, ... , 2^J
//' @details 
//' Used in wvar object.
//' @examples
//' scales_cpp(5)
// [[Rcpp::export]]
arma::vec scales_cpp(unsigned int nb_level){
  // Define scales
  arma::vec scales(nb_level);
  for(unsigned int i=0; i< nb_level;i++){
    scales(i) = pow(2,i+1);
  }
  return scales;
}

//' @title Computes the (MODWT) wavelet covariance matrix
//' @description Calculates the (MODWT) wavelet covariance matrix
//' @param signal_modwt A \code{field<vec>} that contains the modwt decomposition.
//' @param nb_level A \code{integer} that contains the level of decomposition J.
//' @param compute_v A \code{string} that indicates what kind of matrix should be created. Possible options: "diag" or "none"
//' @param robust A \code{boolean} that triggers the use of the robust estimate.
//' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
//' @return A \code{field<mat>} containing the covariance matrix.
//' @examples
//' \dontrun{
//' x=rnorm(100)
//' decomp = modwt(x)
//' V = compute_cov_cpp(decomp$data, decomp$nlevels, compute_v="diag", robust = TRUE, eff=0.6)
//' }
// [[Rcpp::export]]
arma::field<arma::mat> compute_cov_cpp(arma::field<arma::vec> signal_modwt, unsigned int nb_level, std::string compute_v="diag",bool robust = true, double eff=0.6){
    // Compute asymptotic covariance matrix
  arma::field<arma::mat> V(2);
  //arma::mat V = diagmat(arma::ones<arma::vec>(nb_level));
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
      V(0) = diagmat(2 * Aj / signal_modwt(0).n_elem);
      
      //if(robust){
      V(1) = 1/eff*V(0);
      //}
    }
    
    // Compute confidence intervals
    //up_gauss = vmod.col(0) + R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
    //dw_gauss = vmod.col(0) - R::qnorm(1-p, 0.0, 1.0, 1, 0)*sqrt(diagvec(V));
  }
  
  return V;
}



// [[Rcpp::export]]
arma::vec theo_wv(const arma::vec& theta, const std::vector<std::string>& desc,
                                const arma::vec& wv_empir,
                                const arma::vec& tau, int N){
  
  unsigned int num_desc = desc.size();
  unsigned int i_theta = 0;
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
  
  for(unsigned int i = 0; i < num_desc; i++){
    // AR 1
    if(desc[i] == "AR1"){
      double phi = theta(i_theta);
      ++i_theta;
      double sig2 = theta(i_theta);
      ++i_theta;
      
      // Compute theoretical WV
      wv_theo += ar1_to_wv(phi,sig2,tau);
    }
    // RW
    else if(desc[i] == "RW"){
      double sig2 = theta(i_theta);
      ++i_theta;
      wv_theo += rw_to_wv(sig2,tau);
    }
    // DR
    else if(desc[i] == "DR"){
      double drift = theta(i_theta);
      ++i_theta;
      wv_theo += dr_to_wv(drift,tau);
    }
    // WN
    else{
      double sig2 = theta(i_theta);
      ++i_theta;
      wv_theo += wn_to_wv(sig2,tau);
    }
    
  }

  return wv_theo;
}


//////// OPTIM FUNCTIONS


// computes theoretical wv
inline arma::vec theoretical_wv(const arma::vec& theta, const std::vector<std::string>& desc,
                                const arma::vec& wv_empir,
                                const arma::vec& tau, int N){
  
  unsigned int num_desc = desc.size();
  unsigned int i_theta = 0;
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
  
  for(unsigned int i = 0; i < num_desc; i++){
    // AR 1
    if(desc[i] == "AR1"){
      double phi = arma::as_scalar(logit_inv(theta.row(i_theta)));
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



// computes theoretical wv
inline arma::vec theoretical_wv(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                                const arma::vec& wv_empir,
                                const arma::vec& tau, int N){
  
  unsigned int num_desc = desc.size();
  unsigned int i_theta = 0;
  arma::vec wv_theo = arma::zeros<arma::vec>(tau.n_elem);
  
  for(unsigned int i = 0; i < num_desc; i++){
    // AR 1
    if(desc[i] == "AR1"){
      double phi;
      if(model_type == "imu"){
        phi = arma::as_scalar(logit_inv(theta.row(i_theta)));
      }
      else{ // ssm
        phi = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
      }
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
double objFunStarting(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, model_type, wv_empir, tau, N);
  arma::vec standardized = 1-wv_theo/wv_empir;
	// Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

double objFun(const arma::vec& theta, const arma::mat& omega,
                 const std::vector<std::string>& desc, std::string model_type,
                 const arma::vec& wv_empir, const arma::vec& tau, int N){
                   
  arma::vec wv_theo = theoretical_wv(theta, desc, model_type, wv_empir, tau, N);

  // Compute quandratic form
	arma::vec dif = wv_theo - wv_empir;
	return arma::as_scalar(trans(dif)*omega*dif);
}

/// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta, const std::vector<std::string>& desc, std::string model_type,
                          const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        _["desc"] = desc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, const std::vector<std::string>& desc, std::string model_type,
                     const arma::mat& omega, const arma::vec& tau, const arma::vec& wv_empir, int N){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun),
                        _["omega"] = omega,
                        _["desc"] = desc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau,
                        _["N"] = N);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}


//////////// ARMA SECTION


// hiding this function for the moment
double objFunStarting_ARMA(const arma::vec& theta, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
  arma::vec ar;
  arma::vec ma;

  if(p == 0){
    ar = arma::zeros<arma::vec>(0);
  }else{
    ar = theta.rows(0,p-1);
  }
  
  if(q == 0){
    ma = arma::zeros<arma::vec>(0); 
  }else{
    ma = theta.rows(p,p+q-1);
  }
  
  double sigma2 = theta(theta.n_elem-1);
  arma::vec wv_theo = arma_to_wv(ar, ma, tau, sigma2);

  arma::vec standardized = 1-wv_theo/wv_empir;
  // Compute quandratic form
  return arma::as_scalar(trans(standardized)*(standardized));
}

arma::vec Rcpp_OptimStart_ARMA(const arma::vec& theta, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting_ARMA),
                        _["p"] = p,
                        _["q"] = q,
                        _["tau"] = tau,
                        _["wv_empir"] = wv_empir);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

double objFun_ARMA(const arma::vec& theta, const arma::mat& omega,  int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
                   
  arma::vec ar;
  arma::vec ma;
  
  if(p == 0){
    ar = arma::zeros<arma::vec>(0);
  }else{
    ar = pseudo_logit_inv(theta.rows(0,p-1));
  }
  
  if(q == 0){
    ma = arma::zeros<arma::vec>(0); 
  }else{
    ma = pseudo_logit_inv(theta.rows(p,p+q-1));
  }
  
  arma::vec wv_theo = arma_to_wv(ar, ma, tau, exp(theta(theta.n_elem-1)) );
  // Compute quandratic form
  arma::vec dif = wv_theo - wv_empir;
  return arma::as_scalar(trans(dif)*omega*dif);
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim_ARMA(const arma::vec& theta, const arma::mat& omega, int p, int q, const arma::vec& tau, const arma::vec& wv_empir){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun_ARMA),
                        _["omega"] = omega,
                        _["p"] = p,
                        _["q"] = q,
                        _["tau"] = tau,
                        _["wv_empir"] = wv_empir);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}


///// END ARMA

// [[Rcpp::export]]
arma::vec gen_model(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc){
    arma::vec x  = arma::zeros<arma::vec>(N);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
  	for(unsigned int i = 0; i < num_desc; i++){
  	  // AR 1
  	  if(desc[i] == "AR1"){
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
  	  else {
  	    double sigmaWN = theta(i_theta);
  	    ++i_theta;
  	    x += gen_wn(N, sigmaWN);
  	  }
  }  
    
  return x;
}


// [[Rcpp::export]]
arma::vec set_starting_values(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type){
    arma::vec starting  = arma::zeros<arma::vec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
  	  // AR 1
  	  if(desc[i] == "AR1" ){
        if(model_type == "imu"){
  	      starting(i_theta) = arma::as_scalar(logit(theta.row(i_theta)));
        }else{ // ssm
          starting(i_theta) = arma::as_scalar(pseudo_logit(theta.row(i_theta)));
        }
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
arma::rowvec set_result_values(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type){
    arma::rowvec result  = arma::zeros<arma::rowvec>(theta.n_elem);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // AR 1
  	  if(desc[i] == "AR1"){
        if(model_type == "imu"){
          result(i_theta) = arma::as_scalar(logit_inv(theta.row(i_theta)));
        }else{ // ssm
          result(i_theta) = arma::as_scalar(pseudo_logit_inv(theta.row(i_theta)));
        }
  	    ++i_theta;
  	    result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
  	  else {
        result(i_theta) = exp(theta(i_theta));
  	    ++i_theta;
  	  }
  }  
    
  return result;
}


// [[Rcpp::export]]
arma::vec set_starting_values_arma(arma::vec theta, int p, int q){
  theta.rows(0,p+q-1) = pseudo_logit(theta.rows(0,p+q-1));
  theta(theta.n_elem -1) = log(theta(theta.n_elem -1));
  
  return theta;
}


// [[Rcpp::export]]
arma::vec set_result_values_arma(arma::vec theta, int p, int q){    
  theta.rows(0,p+q-1) = pseudo_logit_inv(theta.rows(0,p+q-1));
  theta(theta.n_elem -1) = exp(theta(theta.n_elem -1));
  
  return theta;
}


// [[Rcpp::export]]
arma::vec gmwm_bootstrapper(const arma::vec&  theta, const std::vector<std::string>& desc, 
                            unsigned int tau, unsigned int N, bool robust, double eff,
                            unsigned int B = 100){
  unsigned int nb_level = floor(log2(N));
  	
	arma::mat res(B, tau+1);
	for(unsigned int i=0; i<B; i++){
		arma::vec x = gen_model(N, theta, desc);

    // MODWT transform
    arma::field<arma::vec> signal_modwt = modwt_cpp(x, "haar", nb_level);
    arma::field<arma::vec> signal_modwt_bw = brick_wall(signal_modwt, haar_filter());
  
		arma::vec wv_x = wave_variance(signal_modwt_bw, robust, eff);
    
	  res.row(i) = arma::trans(wv_x);
	}
	return cov(res);
}

inline arma::vec ar1_draw(unsigned int num_ars, double sigma_tot, std::string model_type){
  
  unsigned int num_params = 2*num_ars;
  arma::vec temp(num_params);
  
  for(unsigned int i = 0; i < num_ars; i++){
    
    // Draw from triangle distributions for phi
    double U = R::runif(0.0, 1.0/3.0);
    
    if(i == 0){
      if(model_type == "imu"){
        // Draw for phi
        temp(2*i) = 1.0/5.0*(1.0-sqrt(1.0-3.0*U));
        temp(2*i+1) = R::runif(0.95*sigma_tot*(1-square(temp(2*i))), sigma_tot);
      }
      else{ // ssm
        // Draw for phi
        temp(2*i) = R::runif(-0.9999999999999, 0.9999999999999);
        // Draw for sigma
        temp(2*i+1) = R::runif(0.0000000000001, sigma_tot);
      }
    }
    else{
      
      if(i!=1){
          // Draw for phi on i >= 3
          temp(2*i) = R::runif(temp(2*(i-1)),0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
      }
      else{
          // Draw for phi on i==1
          temp(2*i) = R::runif(0.995,0.9999999); //1.0/40.0*(38.0-sqrt(6.0*U-2.0)) + .05;
      }
      
      // Draw for process variance
      temp(2*i+1) = R::runif(0.0, 0.01*sigma_tot*(1-square(temp(2*i+1))) );
    } // end if
    
  } // end for
  
  return temp;
}

inline arma::vec unif_sigma_sample(unsigned int num, double start, double end){
  arma::vec temp(num);
  
  for(unsigned int i = 0; i<num; i++){
    temp(i) = R::runif(start,end);
  }
  
  return temp;
}


// @title Randomly guess a starting parameter
// @description Sets starting parameters for each of the given parameters. 
// @usage guess_initial(signal, w, desc, model_type, num_param, wv_empir, tau, N, B)
// @param signal A \code{vec} that contains the data
// @param w A \code{map<string,int>} that lists supported models and the amount in the model.
// @param model_type A \code{string} that indicates whether it is an SSM or IMU.
// @param num_params An \code{unsigned int} number of parameters in the model (e.g. # of thetas).
// @param wv_empir A \code{vec} that contains the empirical wavelet variance.
// @param tau A \code{vec} that contains the scales. (e.g. 2^(1:J))
// @param N A \code{integer} that indicates the signal length
// @param B A \code{integer} that indicates how many random draws that should be performed.
// @return A \code{vec} containing smart parameter starting guesses to be iterated over.
// @name guess_initial
// @docType methods
// @rdname guess_initial-methods
arma::vec guess_initial(arma::vec signal, std::map< std::string ,int>& w, const std::vector<std::string>& desc, std::string model_type,
                        unsigned int num_param, const arma::vec& wv_empir, const arma::vec& tau, 
                        unsigned int N, unsigned int B=1000){
  
  // Obtain the sum of variances for sigma^2_total.
  double sigma_tot = arma::sum(wv_empir);
    
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
          temp_theta.rows(i_theta, i_theta + 2*out - 1) = ar1_draw(out, sigma_tot, model_type);
          i_theta += 2*out;
        }
        else if(type == "DR"){   
          double dr_ed = mean(diff_cpp(signal,1));
          if(dr_ed > 0){
            dr_ed = R::runif(0,2*dr_ed);
          }else{
            dr_ed = R::runif(2*dr_ed,0);
          }
          temp_theta.rows(i_theta, i_theta + out - 1).fill( dr_ed  );
          i_theta += out;
        }
        else if(type == "QN"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, .0000001, sigma_tot);
          i_theta += out;
        }
        else if(type == "RW"){
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot/(signal.n_elem*1000.0), 2.0*sigma_tot/signal.n_elem);
          i_theta += out;
        }
        else{ // WN
          temp_theta.rows(i_theta, i_theta + out - 1) = unif_sigma_sample(out, sigma_tot/2.0, sigma_tot);
          i_theta += out;
        }
      } // end if
    } // end for

    double obj = objFunStarting(temp_theta, desc, model_type, wv_empir, tau, N);
    
    if(min_obj_value > obj){
      min_obj_value = obj;
      starting_theta = temp_theta;
    } //end if
  } // end for
  
  return starting_theta;
}

// Counts 
inline std::map<std::string, int> counted_map(const std::vector<std::string>& desc){
  std::map<std::string, int> w;
  w["AR1"]=0;
  w["DR"]=0;
  w["RW"]=0;
  w["QN"]=0;
  w["WN"]=0;

  for (unsigned int i = 0; i < desc.size(); i++) {
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
//' @usage count_models_alpha(desc)
//' @param desc A \code{vector<string>} that contains the model type
//' @return A \code{vec} with the model counts.
//' @details 
//' The types of models supported are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' Note: This order is due to the fact that the models are created by inserting into a C++ Map. The map sorts these values...
//' This function is NOT available in R. 
// [[Rcpp::export]]
arma::vec count_models_alpha(const std::vector<std::string>& desc) {
  
    std::map<std::string, int> w = counted_map(desc);

    arma::vec num_models = arma::zeros<arma::vec>(w.size());
    
    for (std::map<std::string, int>::iterator p = w.begin(); p != w.end(); ++p) {
              int out = p->second;
              num_models(std::distance(w.begin(), p)) = out;
    }
    return num_models;
}

// [[Rcpp::export]]
unsigned int num_model_params(const std::vector<std::string>& desc) {
  // Count the number of models we are working with
  std::map<std::string, int> w = counted_map(desc);
  
  // Return the total number of parameters we need to setup.
  return count_params(w);
}


//' @title User Specified Initial Values for GMWM Estimator
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model.
//' @usage adv_gmwm_imu_cpp(theta, desc, V, wv_empir, tau, N)
//' @param theta A \code{vector} with dimensions N x 1 that contains user-supplied initial values for parameters
//' @param desc A \code{vector<string>} indicating the models that should be considered.
//' @param V A \code{matrix} that represents the covariance matrix.
//' @param wv_empir A \code{vector} that contains the empirical wavelet variance
//' @param N A \code{integer} that indicates the length of the signal being studied.
//' @return A \code{vec} that contains the parameter estimates from GMWM estimator.
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). If type = "IMU" or "SSM", then
//' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
//' options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' If type = "ARMA", the function takes condition least squares as starting values; if type = "IMU" or type = "SSM" then
//' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec adv_gmwm_imu_ssm_cpp(const arma::vec& theta, const std::vector<std::string>& desc, std::string model_type,
                   const arma::mat& V, const arma::vec& wv_empir,
                   const arma::vec& tau, unsigned int N){
                                 
  // Number of parameters
  //unsigned int num_param = theta.n_elem;
      
  // Starting values
  arma::vec starting_theta = set_starting_values(theta, desc, model_type);
  
  // Optimize Starting values via Jannick's Method
  starting_theta = Rcpp_OptimStart(starting_theta, desc, model_type, tau, wv_empir, N);

  // ------------------------------------
  // Compute standard GMWM
  // ------------------------------------
      
  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  // Find GMWM estimator
  arma::vec estim_GMWM = Rcpp_Optim(starting_theta, desc, model_type, omega, tau, wv_empir, N);
  
  return set_result_values(estim_GMWM, desc, model_type);                          
}


//' @title GMWM for IMU and SSM
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model.
//' @usage gmwm_imu_ssm_cpp(desc, signal, model_type, V, wv_empir, tau, N, B = 1000)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param model_type A \code{character string} indicating if the function should estimate an ARMA model ("ARMA"), a model for IMU sensor calibration ("IMU") or a state-space model ("SSM")
//' @param params A \code{vector} being numeric (if type = "ARMA") or character string (if type = "IMU" or type = "SSM")
//' @param robust A \code{bool} indicating if the function should provide a robust estimation of the model parameters (by default = FALSE).
//' @return gmwm A \code{list} that contains:
//' \itemize{
//'  \item{par}{The estimated model parameters}
//'  \item{CI}{The 95\% confidence intervals for the estimated model parameters.}
//' }
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). If type = "IMU" or "SSM", then
//' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
//' options are:
//' \itemize{
//'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
//'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
//'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
//'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
//'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
//' }
//' If type = "ARMA", the function takes condition least squares as starting values; if type = "IMU" or type = "SSM" then
//' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec gmwm_imu_ssm_cpp(const std::vector<std::string>& desc, const arma::vec& signal, std::string model_type,
               const arma::mat& V, const arma::vec& wv_empir,
               const arma::vec& tau, unsigned int N, unsigned int B = 1000){
  
  // Count the number of models we are working with
  std::map<std::string, int> w = counted_map(desc);
  
  // Return the total number of parameters we need to setup.
  unsigned int num_param = count_params(w);
  
  // Give it a guess
  arma::vec guess_me = guess_initial(signal, w, desc, model_type, num_param, wv_empir, tau, N, B);

  // And return value...
  return adv_gmwm_imu_ssm_cpp(guess_me, desc, model_type, V, wv_empir, tau, N);                        
}

//' @title GMWM ARMA
//' @description This function uses the Generalized Method of Wavelet Moments to estimate the parameters of a time series model based on ARMA.
//' @usage gmwm_arma_cpp(theta, V, p, q, tau, wv_empir)
//' @param x A \code{vector} with dimensions N x 1. 
//' @param type A \code{character string} indicating if the function should estimate an ARMA model ("ARMA"), a model for IMU sensor calibration ("IMU") or a state-space model ("SSM")
//' @param params A \code{vector} being numeric (if type = "ARMA") or character string (if type = "IMU" or type = "SSM")
//' @param robust A \code{bool} indicating if the function should provide a robust estimation of the model parameters (by default = FALSE).
//' @return gmwm A \code{list} that contains:
//' \itemize{
//'  \item{par}{The estimated model parameters}
//'  \item{CI}{The 95\% confidence intervals for the estimated model parameters.}
//' }
//' @details
//' The function estimates a variety of time series models. If type = "ARMA" then the parameter vector (param) should
//' indicate the order of the AR process and of the MA process (i.e. param = c(AR,MA)). 
//' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
//' 
//' @author JJB
//' @references Wavelet variance based estimation for composite stochastic processes, S. Guerrier and Robust Inference for Time Series Models: a Wavelet-Based Framework, S. Guerrier
//' @keywords internal
//' @examples
//' # Coming soon
// [[Rcpp::export]]
arma::rowvec gmwm_arma_cpp(const arma::vec& theta, const arma::mat& V, unsigned int p, unsigned int q,
                const arma::vec& tau, const arma::vec& wv_empir){
    
  arma::vec starting_theta = set_starting_values_arma(theta, p, q);

  // Omega matrix
  arma::mat omega = arma::inv(diagmat(V));
  
  arma::vec values = Rcpp_Optim_ARMA(theta, omega, p, q, tau, wv_empir);
  
  // Find GMWM estimator
    // Initialize it
  arma::rowvec GMWM = arma::trans(set_result_values_arma(values, p, q));
      
  return GMWM;                          
}


/*

//' @title Simulate GMWM
//' 
//' @examples
//' x=rnorm(100)
//' wavelet_variance_cpp(x, "haar", "diag")
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
  	arma::field<arma::mat> wv_x = wavelet_variance_cpp(x, "haar", "diag");
  	
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


//
// Mondal and Percival estimator
//

//What/Where is get.slepians?
arma::vec percival(arma::vec x){
  arma::vec xsq = arma::square(x);
  double Tn = log(median(xsq));
  arma::mat beta = get.slepians(npoints=x.n_elem,nwin=5)/x.n_elem;
  arma::rowvec colsums = arma::sum(beta); // 1 x n
  arma::vec J = arma::trans(beta)*arma::sign(log(xsq)-Tn);
  arma::vec mu = (J * colsums)/(colsums*arma::trans(colsums));
  arma::mat Ahat = arma::mean(arma::square(J-mu*colsums));
  
  double temp = R::qnorm(3.0/4.0, 0.0, 1.0, 1,0); //0.6744898
  // dnorm(qnorm(3/4)) = 0.3177766
  
  arma::vec muhat = Tn-2*log(temp)-(Ahat/square(-2*R::dnorm(temp, 0.0, 1.0, 0)*temp))/(2*x.n_elem);
  return exp(muhat);
}
*/