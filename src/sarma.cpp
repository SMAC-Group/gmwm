#include <RcppArmadillo.h>
#include "sarma.h"

arma::vec {
  // Number of ARMA(p,q) parameters
  unsigned int np = objdesc(0), nq = objdesc(1);
  
  // Number of Seasonal (P,Q)
  unsigned int nsp = objdesc(2), nsq = objdesc(3), ns = objdesc(4);
  
  // Find the total number of parameters to expand
  unsigned int p = np + ns * nsp, q = nq + ns * nsq;
}

//' Expand Parameters for an SARMA object
//' 
//' Creates an expanded PHI and THETA vector for use in other objects. 
//' @param params  A \code{vec} containing the theta values of the parameters.
//' @param objdesc A \code{vec} containing the model term information.
//' @return A \code{field<vec>} of size two as follows:
//' \itemize{
//'   \item AR    values
//'   \item THETA values
//' }
//' @details 
//' The \code{objdesc} is assumed to have the structure of:
//' \itemize{
//' \item AR(p)
//' \item MA(q)
//' \item SAR(P)
//' \item SMA(Q)
//' \item Seasons
//' }
//' @keywords internal
//' @examples
//' m = expand_sarima(c(0.5,.2,0,.1,.92,.83,.42,.33,.12), c(2,2,2,3,12))
// [[Rcpp::export]]
arma::field<arma::vec> expand_sarma(const arma::vec& params, const arma::vec& objdesc) {

    
    
    // Loop variables
    unsigned int i, j;
    
    // Fill with zeros in advance.
    arma::vec phi   = arma::zeros<arma::vec>(p);
    arma::vec theta = arma::zeros<arma::vec>(q);
    
    /* Fill the non-seasonal components */
    
    // Fill AR(p)
    for (i = 0; i < np; i++) phi(i) = params(i);
    // Fill MA(q)
    for (i = 0; i < nq; i++) theta(i) = params(i + np);
    
    
    /* Is a seasonal component present? */
    if (ns > 0) {
      
      /* Fill the seasonal components */
      
      // Fill the Seasonal AR(P)
      for (j = 0; j < nsp; j++) {
        phi((j + 1) * ns - 1) += params(j + np + nq);
        
        // Handle the interaction between AR and SAR terms
        for (i = 0; i < np; i++){
          phi((j + 1) * ns + i) += params(i) * params(j + np + nq);
        }
      }
      
      // Fill the Seasonal MA(Q)
      for (j = 0; j < nsq; j++) {
        theta( (j + 1) * ns - 1 ) += params(j + np + nq + nsp);
        
        // Handle interaction between MA and SMA terms
        for (i = 0; i < nq; i++){
            theta( (j + 1) * ns + i ) += params(i + np) * params(j + np + nq + nsp);
        }
      }
    }
    
    // Release transformed parameter space.
    arma::field<arma::vec> out(2);
    
    out(0) = phi;
    out(1) = theta;
    
    return out;
}