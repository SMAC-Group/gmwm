#include <RcppArmadillo.h>

#include "inference.h"

// Need for getObjFun
#include "objective_functions.h"

// Need for gmwm_engine in chisq test
#include "gmwm_logic.h"

using namespace Rcpp;

//' @title Calculate the Psi matrix
//' @description Computes the Psi matrix using supplied parameters
//' @param A first derivative matrix
//' @param v_hat bootstrapped V
//' @param omega original omega matrix
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat calculate_psi_matrix(const arma::mat& A, const arma::mat& v_hat, const arma::mat& omega){ 
  arma::mat A_trans = arma::trans(A);
  arma::mat B = arma::inv(A_trans*omega*A)*A_trans*omega;
  
  return B*v_hat*arma::trans(B);
}

//' @title Format the Confidence Interval for Estimates
//' @description Creates hi and lo confidence based on SE and alpha.
//' @param theta
//' @param se
//' @param alpha
//' @return A \code{mat} that has:
//' \itemize{
//' \item Column 1: Lo CI
//' \item Column 2: Hi CI
//' \item Column 3: SE
//' }
// [[Rcpp::export]]
arma::mat format_ci(const arma::vec& theta,
                    const arma::vec& se,
                    double alpha){
  
  double z = R::qnorm(1-alpha, 0, 1, true, false);
  
  unsigned int nparams = theta.n_elem;
  arma::mat ci(nparams, 3);
  
  ci.col(0) = theta - z*se;
  ci.col(1) = theta + z*se;
  ci.col(2) = se;
  
  return ci;
}

//' @title Generate the Confidence Interval for Theta Estimates
//' @description yaya
//' @param theta
//' @param psi
//' @param z
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat theta_ci(const arma::vec& theta,
                   const arma::mat& A, 
                   const arma::mat& v_hat, const arma::mat& omega, double alpha){

  arma::mat psi = calculate_psi_matrix(A, v_hat, omega);
  
  arma::vec se = sqrt(diagvec(psi));

  return format_ci(theta, se, alpha);
}


//' @title Compute the GOF Test
//' @description yaya
//' @param theta
//' @param psi
//' @param p
//' @return A \code{vec} that has
//' \itemize{
//' \item Test Statistic
//' \item P-Value
//' \item DF
//' } 
// [[Rcpp::export]]
arma::vec gof_test(arma::vec theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   std::string model_type,
                   const arma::vec& tau,
                   const arma::mat& v_hat, const arma::vec& wv_empir){
  
  
  arma::vec estimate = gmwm_engine(theta,
                                   desc, objdesc, 
                                  model_type, 
                                  wv_empir,
                                  v_hat,
                                  tau,
                                  false); // starting is false

  double test_stat = getObjFun(estimate, desc, objdesc, model_type,
                                arma::inv(v_hat), wv_empir, tau);
    
  unsigned int df = tau.n_elem - theta.n_elem;
  
  double p_value = 1.0 - R::pchisq(test_stat, df, true, false);

  arma::vec out(3);
  out(0) = test_stat;
  out(1) = p_value;
  out(2) = df;

  return out;
}
