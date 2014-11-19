#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Compute Allan Variance
//' @description Computation of Allan Variance
//' @usage allan_variance(x)
//' @param x A \code{vector} with dimensions N x 1. 
//' @return av A \code{matrix} that contains the cluster size as the first column and the allan variance as the second.
//' @author JJB
//' @examples
//' white.noise = rnorm(500, 0, 1)
//' allan_variance(white.noise)
// [[Rcpp::export]]
arma::mat allan_variance(arma::vec x) {
  
   // Length of vector
   unsigned int T = x.n_elem;
   // Create the number of halves possible and use it to find the number of clusters
   double J = floor(log10( ceil((T-1)/2) )/log10(2));
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J+1,2);
   
   for (int i = 0; i <= J; i++){
     // Tau / cluster size.
     unsigned int tau = pow(2,i);
     
     // Obtain the number of the clusters
     unsigned int nclusters = floor(T/tau);
     arma::vec averages = arma::zeros<arma::vec>( nclusters ); 
     
     // Get the averages
     unsigned int row_index = 0;
     for(unsigned int cluster_num = 0; cluster_num < nclusters; cluster_num++){
       /* Benchmark difference between sum and accu? */
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
  }
  
  return av;
}
