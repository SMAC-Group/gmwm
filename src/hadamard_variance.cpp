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


#include "hadamard_variance.h"


/* ----------------------- Start Hadamard Variance Functions ------------------------ */

//' @title Compute Tau-Overlap Hadamard Variance
//' @description Computation of  Hadamard Variance
//' @usage hadam_to_cpp(x)
//' @param x A \code{vector} with dimensions M x 1. 
//' @return  A \code{matrix} that contains:
//' \itemize{
//'  \item{Col 1}{The size of the cluster}
//'  \item{Col 2}{The Hadamard variance}
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
//' @keywords internal
// [[Rcpp::export]]
arma::mat hadam_to_cpp(arma::vec x) {
  
  // Length of vector
  unsigned int T = x.n_elem;
  // Create the number of halves possible and use it to find the number of clusters
  unsigned int J = floor(log10(T)/log10(3))-2;
  
  // Hadamard Variance Matrix
  arma::mat hv = arma::zeros<arma::mat>(J,3);
  
  for (unsigned int i = 1; i <= J; i++){
    // Tau
    unsigned int tau = pow(3,i);
    
    // Y.Bar
    unsigned int N = floor(T/tau);
    arma::vec yBar = arma::zeros<arma::vec>(N);
    for(unsigned int j = 0; j < N;j++){
      yBar(j) = sum( x.rows(tau*j, tau*j+tau - 1) )/tau;
    }

    // Clusters
    unsigned int M = floor(T/(3*tau)) ;
    

    double summed = 0;
    for(unsigned int k = 0; k <=  M - 1 ; k++){
      summed +=  pow( yBar(3*k+2) - 2*yBar(3*k+1) + yBar(3*k),2);
    }

    // Cluster size
    hv(i-1,0) = tau; 
    // Compute the Hadamard Variance estimate
    hv(i-1,1) = summed/(6*(M-2)); 
    // Compute Error
    hv(i-1,2) = 1/sqrt(2*( (double(T)/tau) - 1) );
  }
  
  // Return as list
  return hv;
}

//' @title Compute Maximal-Overlap Hadamard Variance using Means
//' @description Computation of Maximal-Overlap Hadamard Variance
//' @usage avar_mo_cpp(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"hadamard"}{The Hadamard variance}
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
//' @keywords internal
// [[Rcpp::export]]
arma::mat hadam_mo_cpp(arma::vec x) {
  
  // Length of vector
  unsigned int T = x.n_elem;
  
  // Create the number of halves possible and use it to find the number of clusters
  unsigned int J = floor(log10(T)/log10(3))-1;
  
  // Hadamard Variance Matrix
  arma::mat ha = arma::zeros<arma::mat>(J,3);
  
  for (unsigned int i = 1; i <= J; i++){
    // Tau
    unsigned int tau = pow(3,i);
    
    // Y.Bar
    arma::vec yBar = arma::zeros<arma::vec>(T);

    for(unsigned int j = 0; j <= T - tau; j++){
      yBar(j) = sum( x.rows(j, j+tau-1) ) / tau;
    }
    
    // Clusters
    unsigned int M = T-3*tau;
    double summed = 0;

    for(unsigned int k = 0; k <= M - 1; k++){
      summed +=  pow( yBar(k+2*tau) - 2*yBar(k+tau) + yBar(k),2);
    }
    
    // Cluster size
    ha(i-1,0) = tau; 
    // Compute the Hadamard Variance estimate
    ha(i-1,1) = summed/(6*(T - 3*tau +1)); 
    // Compute Error
    ha(i-1,2) = 1/sqrt(2*( (double(T)/tau) - 1) );
  }
  
  return ha;
}

/* --------------------- End Hadamard Variance Functions ---------------------- */