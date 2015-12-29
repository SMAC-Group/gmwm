
/* Copyright (C) 2014 - 2015  James Balamuta
 *
 * This file is part of GMWM R Methods Package
 *
 * The file uses methods in the r-to-armadillo project and is free software: you can redistribute it and/or modify it
 * under the terms of the MIT License.
 *
 * The r-to-armadillo project is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */


#include <RcppArmadillo.h>

#include "rtoarmadillo.h"
#include "armadillo_manipulations.h"

/* ----------------- R to Armadillo Functions ------------------ */

// A special define is included in rtoarmadillo.h used in these functions...




//' @title Generate a sequence of values
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param a An \code{int}, that denotes the starting point.
//' @param b An \code{int}, that denotes the ending point.
//' @return A \code{vector} containing values moving from a to b. There are no restrictions on A's range.
//' @author James J Balamuta
//' @keywords internal
//' @examples
//' #Call with the following data:
//' seq_cpp(3, 5)
//' seq_cpp(5, 3)
// [[Rcpp::export]]
arma::vec seq_cpp(int a, int b){
  int d = abs(b-a)+1;
  
  int inc = ( a < b ? 1 : -1 );
  arma::vec s(d);
  
  s.fill(inc);
  s(0) = a;

  return cumsum(s);
}

//' @title Generate a sequence of values based on supplied number
//' @description Creates a vector containing a sequence of values starting at 1 and going to the terminal point.
//' @param n An \code{int} that denotes the length of the vector.
//' @return A \code{vector} containing values moving from 1 to n.
//' @author James J Balamuta
//' @keywords internal
//' @examples 
//' #Call with the following data:
//' seq_len_cpp(5)
// [[Rcpp::export]]
arma::vec seq_len_cpp(unsigned int n){
  arma::vec seq = arma::ones<arma::vec>(n);
  return cumsum(seq);
}


//' @title Find Quantiles
//' @description Attempts to find quantiles
//' @param x A \code{vec} of data
//' @param probs A \code{vec} of the quantiles to find.
//' @return A \code{vector} containing the quantiles
//' @author James J Balamuta
//' @keywords internal
//' @examples 
//' #Call with the following data:
//' quantile_cpp(c(1,2,3,4,5,6,7), c(.25,.5,.75))
//' quantile(c(1,2,3,4,5,6,7), c(.25,.5,.75))
// [[Rcpp::export]]
arma::vec quantile_cpp(arma::vec x, const arma::vec& probs) {
  
  unsigned int n = x.n_elem;
  
  arma::uvec index = arma::conv_to<arma::uvec>::from( (n - 1) * probs);
  arma::uvec lo = floor(index);
  arma::uvec hi = ceil(index);
  
  // bad for sorting large data. need partial sort.
  x = sort(x);
  
  arma::vec qs = x(lo);
  arma::uvec i = index > lo;
  arma::uvec h = (index - lo);
  h = h.elem(i);
  qs.elem(i) = (1 - h) % qs.elem(i) + h % x.elem(hi.elem(i));
  
  return qs;
}


//' @title Lagged Differences in Armadillo
//' @description Returns the ith difference of a time series of rth lag.
//' @param x A \code{vec} that is the time series
//' @param lag A \code{unsigned int} that indicates the lag
//' @param differences A \code{dif} that indicates how many differences should be taken
//' @return A \code{vector} containing the differenced time series.
//' @author JJB
//' @keywords internal
//' @examples
//' x = rnorm(10000, 0, 1)
//' diff_cpp(x,1,1)
// [[Rcpp::export]]
arma::vec diff_cpp(arma::vec x, unsigned int lag, unsigned int differences){
  
  // Difference the series i times
  for(unsigned int i=0; i < differences; i++){
    // Each difference will shorten series length
    unsigned int n=x.n_elem;
    // Take the difference based on number of lags
    x = (x.rows(lag,n-1) - x.rows(0,n-lag-1));
  }
  
  // Return differenced series:
  return x;
}

//' @title Converting an ARMA Process to an Infinite MA Process
//' @description Takes an ARMA function and converts it to an infinite MA process.
//' @param ar A \code{column vector} of length p
//' @param ma A \code{column vector} of length q
//' @param lag_max A \code{int} of the largest MA(Inf) coefficient required.
//' @return A \code{column vector} containing coefficients
//' @details This function is a port of the base stats package's ARMAtoMA. There is no significant speed difference between the two.
//' @author R Core Team and JJB
//' @keywords internal
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
//' @param x A \code{column vector} of length T
//' @param filter A \code{column vector} of length f
//' @param sides An \code{int} that takes either 1:for using past values only or 2: filter coefficients are centered around lag 0.
//' @param circular A \code{bool} that indicates if the filter should be wrapped around the ends of the time series.
//' @return A \code{column vec} that contains the results of the filtering process.
//' @details This is a port of the cfilter function harnessed by the filter function in stats. 
//' It is about 5-7 times faster than R's base function. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
//' @author R Core Team and JJB
//' @keywords internal
//' @examples
//' x = 1:100
//' # 
//' cfilter(x, rep(1, 3), sides = 2, circular = FALSE)
//' # Using R's function
//' filter(x, rep(1, 3))
//' #
//' cfilter(x, rep(1, 3), sides = 1, circular = FALSE)
//' # Using R's function
//' filter(x, rep(1, 3), sides = 1)
//' #
//' cfilter(x, rep(1, 3), sides = 1, circular = TRUE)
//' # Using R's function
//' filter(x, rep(1, 3), sides = 1, circular = TRUE)
// [[Rcpp::export]]
arma::vec cfilter(arma::vec x, arma::vec filter, int sides, bool circular)
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
//' @keywords internal
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

// @title Expand Grid for Same Dimensional Case
// @description Creates the different pairings possible with two different variables.
// @usage expand_grid_red(nx)
// @param nx An \code{integer} of length f that contains the initial values of the time series in reverse.
// @return x A \code{matrix} listing values from 1...nx in one column and 1...1, 2...2,....,n...n, in the other
// @author JJB
// @details This function is hidden and is not accessible from R.
// @name expand_grid_red
// @docType methods
// @rdname expand_grid_red-methods
// @keywords internal
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
//' @keywords internal
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


//' @title Discrete Fourier Transformation for Autocovariance Function
//' @description Calculates the autovariance function (ACF) using Discrete Fourier Transformation.
//' @param x A \code{cx_vec}. 
//' @return A \code{vec} containing the ACF.
//' @details 
//' This implementation is 2x as slow as Rs. 
//' Two issues: 1. memory resize and 2. unoptimized fft algorithm in arma.
//' Consider piping back into R and rewrapping the object. (Decrease of about 10 microseconds.)
//' @keywords internal
//' @examples
//' x=rnorm(100)
//' dft_acf(x)
// [[Rcpp::export]]
arma::vec dft_acf(const arma::vec& x){
    unsigned int n = x.n_elem;
    arma::cx_vec ff = arma::fft(x,2*n);
    
    arma::cx_vec iff = arma::conv_to< arma::cx_vec >::from( arma::square(arma::real(ff)) + arma::square(arma::imag(ff)) ); //expensive
    
    iff = arma::ifft(iff);
  
    iff = iff.rows(0,n-1) / n;
    
    return arma::real(iff);
}

//' @title Mean of the First Difference of the Data
//' @description The mean of the first difference of the data
//' @param x A \code{vec} containing the data 
//' @return A \code{double} that contains the mean of the first difference of the data.
//' @keywords internal
//' @examples
//' x=rnorm(100)
//' mean_diff(x)
// [[Rcpp::export]]
double mean_diff(const arma::vec& x){
  return arma::mean(diff_cpp(x, 1, 1));
}



/* ------------------ End R to Armadillo Functions ----------------------- */