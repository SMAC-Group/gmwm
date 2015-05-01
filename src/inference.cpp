#include <RcppArmadillo.h>

// Need for getObjFun
#include "objective_functions.h"

// Need for gmwm_engine in chisq test
#include "gmwm_logic.h"


//#include "inference.h"
using namespace Rcpp;

//' @title Calculate the Psi matrix
//' @description Computes the Psi matrix using supplied parameters
//' @param D
//' @param v_hat
//' @param omega
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat calculate_psi_matrix(const arma::mat& D, const arma::mat& v_hat, const arma::mat& omega){ 
  arma::mat D_trans = arma::trans(D);
  arma::mat B = arma::inv(D_trans*omega*D)*D_trans*omega;
  
  return B*v_hat*arma::trans(B);
}

//' @title Generate the Confidence Interval for Theta Estimates
//' @description yaya
//' @param theta
//' @param psi
//' @param p
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::mat theta_ci(arma::vec theta, arma::mat psi, double alpha){
  unsigned int nparams = theta.n_elem;
  arma::mat ci(nparams, 3);
  double z = R::qnorm(1-alpha, 0, 1, true, false);
  
  arma::vec se = sqrt(diagvec(psi));
  
  ci.col(0) = theta - z*se;
  ci.col(1) = theta + z*se;
  ci.col(2) = se;

  return ci;
}

//' @title Compute the GOF Test
//' @description yaya
//' @param theta
//' @param psi
//' @param p
//' @return A \code{mat} that has the first column 
// [[Rcpp::export]]
arma::vec gof_test(const arma::vec& theta, 
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
  
  double p_value = R::pchisq(test_stat, df, true, false);

  arma::vec out(3);
  out(0) = test_stat;
  out(1) = p_value;
  out(2) = df;

  return out;
}

// [[Rcpp::export]]
arma::field<arma::mat> inference_summary(const arma::vec& theta, 
                                        const std::vector<std::string>& desc,
                                        const arma::field<arma::vec>& objdesc,
                                        std::string model_type,
                                        const arma::vec& tau,
                                        arma::mat D, 
                                        arma::mat v_hat, arma::mat omega, arma::vec wv_empir, double alpha){
  
  
  
  arma::mat psi = calculate_psi_matrix(D, v_hat, omega);
  arma::mat ci = theta_ci(theta, psi, alpha);
  
  arma::vec gof = gof_test(theta, desc, objdesc, model_type, tau, v_hat, wv_empir);

  arma::field<arma::mat> out(2);
  
  out(0) = ci;
  out(1) = gof;

  return out;
} 