#include <RcppArmadillo.h>

#include "gen_process.h"

// Need to have access to diff_cpp
#include "rtoarmadillo.h"

using namespace Rcpp;
/* ------------------------------ Start Process Generation Functions ------------------------------ */

//' @title Generate a white noise process
//' @description Generates a white noise process with variance parameter sigma.
//' @param N An \code{integer} for signal length.
//' @param sigma2 A \code{double} that contains process variance.
//' @return wn A \code{vec} containing the white noise.
//' @examples
//' gen_wn(10, 1.5)
// [[Rcpp::export]]
arma::vec gen_wn(const unsigned int N, const double sigma2 = 1)
{
  arma::vec wn(N);
  double sigma = sqrt(sigma2);
  for(unsigned int i = 0; i < N; i++){
      wn(i) = R::rnorm(0.0, sigma);
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
arma::vec gen_dr(const unsigned int N, const double slope = 5)
{
  arma::vec gd(N);
  gd.fill(slope);
	return cumsum(gd);
}

//' @title Generate a Quantisation Noise (QN) sequence
//' @description Generate an QN sequence given q2
//' @param N An \code{integer} for signal length.
//' @param q2 A \code{double} that contains autocorrection.
//' @return  A \code{vec} containing the QN process.
//' @details 
//' To generate the quantisation noise, we follow this recipe:
//' First, we generate using a random uniform distribution:
//' \deqn{U_k^*\sim U\left[ {0,1} \right]}{U_k^*~U[0,1]}
//' 
//' Then, we multiple the sequence by \eqn{\sqrt{12}}{sqrt(12)} so:
//' \deqn{{U_k} = \sqrt{12} U_k^*}{U_k = sqrt(12)*U_k^*}
//' 
//' Next, we find the derivative of \eqn{{U_k}}{U_k}
//' \deqn{{{\dot U}_k} = \frac{{{U_{k + \Delta t}} - {U_k}}}{{\Delta t}}}{U_k^. = (U_(k + (delta)t) - U_k)}
//'
//' In this case, we modify the derivative such that:
//' \eqn{{{\dot U}_k}\Delta t = {U_{k + \Delta t}} - {U_k}}{U_k^. * (delta)t = U_{k + (delta)*t} - U_k}
//'
//' Thus, we end up with:
//' \deqn{{x_k} = \sqrt Q {{\dot U}_k}\Delta t}{x_k = sqrt(Q)*U_k^.*(delta)t}
//' \deqn{{x_k} = \sqrt Q \left( {{U_{k + 1}} - {U_k}} \right)}{x_k = sqrt(Q)* (U_(k+1) - U_(k))}
//'
//' @examples
//' gen_qn(10, 5)
// [[Rcpp::export]]
arma::vec gen_qn(const unsigned int N, double q2 = .1)
{
  double sqrt12 = sqrt(12);
  
  arma::vec gu(N+1);
  
	for(unsigned int i=0; i <= N; i++ )
	{		
		gu(i) = sqrt12*R::runif(0.0,1.0);
	}

	return sqrt(q2)*diff_cpp(gu);
}


//' @title Generate an AR(1) sequence
//' @description Generate an AR sequence given phi and sig2.
//' @details This needs to be extended to AR(p) see \code{arima.sim} and \code{filter}.
//' @param N An \code{integer} for signal length.
//' @param phi A \code{double} that contains autocorrection.
//' @param sigma2 A \code{double} that contains process variance.
//' @return gm A \code{vec} containing the AR(1) process.
//' @examples
//' gen_ar1(10, 5, 1.2)
// [[Rcpp::export]]
arma::vec gen_ar1(const unsigned int N, const double phi = .3, const double sigma2 = 1)
{

	arma::vec wn = gen_wn(N+1, sigma2);
	arma::vec gm = arma::zeros<arma::vec>(N+1);
	for(unsigned int i=1; i <= N; i++ )
	{		
		gm(i) = phi*gm(i-1) + wn(i);
	}

	return gm.rows(1,N);
}

//' @title Generate a random walk without drift
//' @description Generates a random walk without drift.
//' @param N An \code{integer} for signal length.
//' @param sigma2 A \code{double} that contains process variance.
//' @return grw A \code{vec} containing the random walk without drift.
//' @examples
//' gen_rw(10, 8.2)
// [[Rcpp::export]]
arma::vec gen_rw(const unsigned int N, const double sigma2 = 1)
{
  arma::vec grw(N);
  double sigma = sqrt(sigma2);
  for(unsigned int i = 0; i < N; i++){
      grw(i) = R::rnorm(0.0, sigma);
  }
  return cumsum(grw);
}


/// [[Rcpp::export]]
arma::vec gen_model(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc){
    arma::vec x  = arma::zeros<arma::vec>(N);
    unsigned int i_theta = 0;
    unsigned int num_desc = desc.size();
    
    for(unsigned int i = 0; i < num_desc; i++){
      // Need to add ARMA generation
      
      double theta_value = theta(i_theta);
  	  // AR 1
  	  if(desc[i] == "AR1"){
  	    ++i_theta;
  	    double sig2 = theta(i_theta);
  	    
  	    // Compute theoretical WV
  	    x += gen_ar1(N, theta_value, sig2);
  	  }
      else if(desc[i] == "ARMA"){
        // Add! via arima.sim
      }
      // DR
  	  else if(desc[i] == "DR"){
  	    x += gen_dr(N, theta_value);
  	  }
      // QN
  	  else if(desc[i] == "QN"){
  	    x += gen_qn(N, theta_value);
  	  }
      // RW
  	  else if(desc[i] == "RW"){
  	    x += gen_rw(N, theta_value);
  	  }
  	  // WN
  	  else {
  	    x += gen_wn(N, theta_value);
  	  }
      
      ++i_theta;
  }  
    
  return x;
}



/* --------------------- END Process Generation Functions -------------------------- */
