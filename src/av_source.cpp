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
//' @references 
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
   unsigned int J = floor(log10( ceil(T/2) )/log10(2));
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 0; i < J; i++){
     // n / cluster size.
     unsigned int n = pow(2,i);
     
     // Obtain the number of the clusters
     unsigned int nclusters = floor((T-1)/n) - 1; 
     arma::vec averages = arma::zeros<arma::vec>( nclusters - 1); 
     
     // Get the averages
     unsigned int row_index = 0;
     for(unsigned int cluster_num = 0; cluster_num < nclusters - 1; cluster_num++){
       averages(cluster_num) = sum(x.rows( row_index, row_index + n - 1)) /n;
       row_index = row_index + n;
     }
     
     Rcpp::Rcout << "Iteration: " << i << std::endl;
     Rcpp::Rcout << "Clusters: " << nclusters << std::endl;

                    
     // Get the difference of the averages
     double summed = 0;
     for (unsigned int k = 0; k <  nclusters - 2; k++){
       summed += pow( (averages(k+1)-averages(k)), 2);
     }
     
     // Cluster size
     av(i,0) = n; 
     // Compute the Allan Variance estimate
     av(i,1) = summed/(2*(nclusters));
     // Compute Error
     av(i,2) = 1/sqrt(2*( (double(T)/n) - 1) );
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
//' \deqn{{Z^2}\left( {i\tau ,\tau } \right) = {\left( {{X_{\tau i}} - 2{X_{\tau i + \tau }} + {X_{\tau i + 2\tau }}} \right)^2}}{See the PDF Documentation}
//' where \eqn{{Z^2}\left( {t\tau ,\tau } \right) = {\left( {{X_t} - 2{X_{t + \tau }} + {X_{t + 2\tau }}} \right)^2}}
//' 
//' @author JJB
//' @references Recipes for Degrees of Freedom of Frequency Stability Estimators, Charles A. Greenhall
//' @references 
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avari_to_arma(combined.ts)
// [[Rcpp::export]]
Rcpp::List avari_to_arma(arma::vec x) {
   // Length of vector
   unsigned int T = x.n_elem;
   // Create the number of halves possible and use it to find the number of clusters
   unsigned int J = floor(log10( ceil(T/2) )/log10(2));
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 0; i < J; i++){
     // n / cluster size.
     unsigned int n = pow(2,i);
     // Saved
     unsigned int n2 = 2*n; // 2 * 2^i
     
     // Obtain the number of the clusters
     unsigned int nclusters = floor((T-1)/n) - 1; 
    
     Rcpp::Rcout << "Iteration: " << i << std::endl;
     Rcpp::Rcout << "Clusters: " << nclusters << std::endl;

     double summed = 0;
     for (unsigned int k = 0; k <  nclusters - 1; k++){
       summed += pow( (  x(k*n+n2) - 2*x(k*n+n) + x(k*n)), 2);
     }
     
     // Cluster size
     av(i,0) = n; 
     // Compute the Allan Variance estimate
     av(i,1) = summed/(2*pow(n,2)*(nclusters - 1));
     // Compute Error
     av(i,2) = 1/sqrt(2*( (double(T)/n) - 1) );
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


//' @title Compute Maximal-Overlap Allan Variance using Iterative Method
//' @description Computation of Maximal-Overlap Allan Variance
//' @usage allan_variance(x)
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
//' \deqn{\sigma _y^2\left( \tau  \right) = V\left( {M,\tau ,{\tau _0}} \right) = \frac{1}{{2{\tau ^2}M}}\sum\limits_{i = 0}^{M - 1} {{Z^2}\left( {i\tau_0 ,\tau } \right)}}
//' where \eqn{{Z^2}\left( {t\tau_0 ,\tau } \right) = {\left( {{X_t} - 2{X_{t + \tau }} + {X_{t + 2\tau }}} \right)^2}}, note \eqn{\tau_0 = 1}.
//' 
//' @author JJB
//' @references Recipes for Degrees of Freedom of Frequency Stability Estimators, Charles A. Greenhall
//' @references 
//' @examples
//' set.seed(999)
//' # Simulate white noise (P 1) with sigma^2 = 4
//' N = 100000
//' white.noise = rnorm(N, 0, 2)
//' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
//' #Simulate random walk (P 4)
//' random.walk = cumsum(0.1*rnorm(N, 0, 2))
//' combined.ts = white.noise+random.walk
//' av_mat = avari_mo_arma(combined.ts)
// [[Rcpp::export]]
Rcpp::List avari_mo_arma(arma::vec x) {
  
   // Length of vector
   unsigned int T = x.n_elem;
   // Create the number of halves possible and use it to find the number of clusters
   unsigned int J = floor(log10( T )/log10(2));
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 0; i < J; i++){
     // Tau / cluster size.
     unsigned int n = pow(2,i);
     // Saved
     unsigned int n2 = 2*n;
     
     // Obtain the number of the clusters
     unsigned int nclusters = T-n2; 
     
     double summed = 0;
     for (unsigned int k = 0; k <  nclusters - 1; k++){
       summed += pow( (x(k) - 2*x(k+n) + x(k+n2) ), 2);
     }
     
     // Cluster size
     av(i,0) = n; 
     // Compute the Allan Variance estimate
     av(i,1) = summed/(2*pow(n,2)*(nclusters));
     // Compute Error
     av(i,2) = 1/sqrt(2*( (double(T)/n) - 1) );
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



//' @title Compute Maximal-Overlap Allan Variance using Means
//' @description Computation of Maximal-Overlap Allan Variance
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
//' \deqn{\frac{1}{{2\left( {N - 2\tau  + 1} \right)}}\sum\limits_{i = 2\tau }^N {\left( {{{\bar y}_i}\left( \tau  \right) - {{\bar y}_{i - \tau }}\left( \tau  \right)} \right)}}{See the PDF documentation}
//' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
//' 
//' 
//' @author JJB
//' @references Recipes for Degrees of Freedom of Frequency Stability Estimators, Charles A. Greenhall
//' @references 
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
   unsigned int J = floor(log10( ceil(T/2) )/log10(2));
   
   // Allan Variance Matrix
   arma::mat av = arma::zeros<arma::mat>(J,3);
   
   for (unsigned int i = 0; i < J; i++){
     // n / cluster size.
     unsigned int n = pow(2,i);
     
     unsigned int n2 = 2*n;
     
     // Obtain the number of the clusters
     unsigned int nclusters = T-n2;
     arma::vec averages = arma::zeros<arma::vec>( nclusters - 1); 
     
     // Get the averages
     unsigned int row_index = 0;
     for(unsigned int cluster_num = 0; cluster_num < nclusters - 1; cluster_num++){
       averages(cluster_num) = sum(x.rows( row_index, row_index + n - 1)) /n;
       row_index = row_index + 1; // increment by 1
     }
     
     Rcpp::Rcout << "Iteration: " << i << std::endl;
     Rcpp::Rcout << "Clusters: " << nclusters << std::endl;

                    
     // Get the difference of the averages
     double summed = 0;
     for (unsigned int k = 0; k <  nclusters - 2; k++){
       summed += pow( (averages(k+1)-averages(k)), 2);
     }
     
     // Cluster size
     av(i,0) = n; 
     // Compute the Allan Variance estimate
     av(i,1) = summed/(2*(nclusters + 1)); // Add one
     // Compute Error
     av(i,2) = 1/sqrt(2*( (double(T)/n) - 1) );
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
