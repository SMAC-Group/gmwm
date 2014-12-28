#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Compute Tau-Overlap Allan Variance
//' @description Computation of Tau-Overlap Allan Variance
//' @usage avarto_arma(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"allan"}{The allan variance}
//'  \item{"errors"}{The error associated with the variance calculation.}
//' }
//' @details
//' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
//' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
//' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
//' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist. 
//' The tau-overlap estimator is given as:
//' 
//' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
//' 
//' @author JJB
//' @references Recipes for Degrees of Freedom of Frequency Stability Estimators, Charles A. Greenhall
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avar_to_arma(combined.ts)
// [[Rcpp::export]]
Rcpp::List avar_to_arma(arma::vec x) {
  
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
  
  
  // Prep export
  arma::vec clusters = av.col(0);
  arma::vec allan = av.col(1);
  arma::vec errors = av.col(2);
  
  // Return as list
  return Rcpp::List::create(
            Rcpp::Named("clusters", clusters),
            Rcpp::Named("allan", allan),
            Rcpp::Named("errors", errors)
          );
}

/*
# Nonoverlapped estimator given in (4)
avar = function(y){
  N = length(y)
	J = floor(log2(length(y))) - 1
	av = rep(NA,J)
	for (i in 1:J){
		tau = 2^i
		yBar = y.bar(y,tau)
		M = floor( N/(2*tau) )
    print(paste("Value of M:", M))
		summed = rep(NA,M)
		for (k in 1:M){
			summed[k] =  (yBar[(2*k)] - yBar[(2*k - 1)])^2
		}
		av[i] = sum(summed)/(2*M)
	}
	return(av)
}

y.bar = function(y,tau){
  N = floor(length(y)/tau)
  yBar = rep(NA,N)
  for (i in 1:N){
    index = (1:tau)+(i-1)*tau
    yBar[i] = mean(y[index])
  }
  return(yBar)
}


# Maximal-overlap estimator given in (5)
avar2 = function(y){
	N = length(y)
	J = floor(log2(N)) - 1
	av = rep(NA,J)
	for (i in 1:J){
		tau = 2^i
		yBar = y.bar2(y,tau)
		M = (N-2*tau)
		summed = rep(NA,M)
    print(paste("Value of M:", M))
		for (k in 1:M){
			summed[k] =  (yBar[k] - yBar[(k+tau)])^2
		}
		av[i] = sum(summed)/(2*(N - 2*tau + 1))
	}
	return(av)
}

y.bar2 = function(y,tau){
  N = length(y)
  yBar = rep(NA,N)
  for (i in 1:N){
    index = (1:tau)+(i-1)
    yBar[i] = mean(y[index])
  }
  return(yBar)
}
*/

//' @title Compute Maximal-Overlap Allan Variance using Means
//' @description Computation of Maximal-Overlap Allan Varianc e
//' @usage avar_mo_arma(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"allan"}{The allan variance}
//'  \item{"errors"}{The error associated with the variance calculation.}
//' }
//' @details
//' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
//' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
//' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
//' Then, \eqn{M = N - 2n} samples exist. 
//' The Maximal-overlap estimator is given as:
//' \eqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }
//' 
//' @author JJB
//' @references Recipes for Degrees of Freedom of Frequency Stability Estimators, Charles A. Greenhall
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avar_to_arma(combined.ts)
// [[Rcpp::export]]
Rcpp::List avar_mo_arma(arma::vec x) {
  
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
  
  
  // Prep export
  arma::vec clusters = av.col(0);
  arma::vec allan = av.col(1);
  arma::vec errors = av.col(2);
  
  // Return as list
  return Rcpp::List::create(
            Rcpp::Named("clusters", clusters),
            Rcpp::Named("allan", allan),
            Rcpp::Named("errors", errors)
          );
}