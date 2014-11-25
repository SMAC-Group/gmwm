#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Compute Allan Variance
//' @description Computation of Allan Variance
//' @usage allan_variance(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{list} that contains:
//' \itemize{
//'  \item{"clusters"}{The size of the cluster}
//'  \item{"allan"}{The allan variance}
//'  \item{"errors"}{The error associated with the variance calculation.}
//' }
//' @author JJB
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = allan_variance(combined.ts)
// [[Rcpp::export]]
Rcpp::List avar_arma(arma::vec x) {
  
   // Length of vector
   unsigned int T = x.n_elem;
   // Create the number of halves possible and use it to find the number of clusters
   unsigned int J = floor(log10( ceil((T-1)/2) )/log10(2)) + 1;
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 0; i < J; i++){
     // Tau / cluster size.
     unsigned int tau = pow(2,i);
     
     // Obtain the number of the clusters
     unsigned int nclusters = floor(T/tau);
     arma::vec averages = arma::zeros<arma::vec>( nclusters ); 
     
     // Get the averages
     unsigned int row_index = 0;
     for(unsigned int cluster_num = 0; cluster_num < nclusters; cluster_num++){
       averages(cluster_num) = sum(x.rows( row_index, row_index + tau - 1))  / tau ;
       row_index = row_index + tau;
     }
                    
     // Get the difference of the averages
     double summed = 0;
     for (unsigned int k = 0; k <  nclusters - 1; k++){
       summed += pow( (averages(k+1)-averages(k)), 2);
     }
     
     // Cluster size
     av(i,0) = tau; 
     // Compute the Allan Variance estimate
     av(i,1) = summed/(2*(nclusters - 1));
     // Compute Error
     av(i,2) = 1/sqrt(2*( (double(T)/tau) - 1) );
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