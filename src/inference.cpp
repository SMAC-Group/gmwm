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

#include "inference.h"

#include "bootstrappers.h"

// Need for getObjFun
#include "objective_functions.h"

// Need for gmwm_engine in chisq test
#include "gmwm_logic.h"

//' @title Calculate the Psi matrix
//' @description Computes the Psi matrix using supplied parameters
//' @param A first derivative matrix
//' @param v_hat bootstrapped V
//' @param omega original omega matrix
//' @backref src/inference.cpp
//' @backref src/inference.h
//' @return A \code{mat} that has the first column 
//' @keywords internal
// [[Rcpp::export]]
arma::mat calculate_psi_matrix(const arma::mat& A, const arma::mat& v_hat, const arma::mat& omega){ 
  arma::mat A_trans = arma::trans(A);
  arma::mat B = arma::pinv(A_trans*omega*A)*A_trans*omega;
  
  return B*v_hat*arma::trans(B);
}

//' @title Format the Confidence Interval for Estimates
//' @description Creates hi and lo confidence based on SE and alpha.
//' @param theta A \code{vec} containing the estimates
//' @param se A \code{vec} containing the standard error
//' @param alpha A \code{double} that contains the confidence level.
//' @return A \code{mat} that has:
//' \itemize{
//' \item Column 1: Lo CI
//' \item Column 2: Hi CI
//' \item Column 3: SE
//' }
//' @backref src/inference.cpp
//' @backref src/inference.h
//' @keywords internal
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
//' @description Create an Asymptotic CI for the Theta Estimates.
//' @param theta A \code{vec} containing the estimates
//' @param A A \code{mat} that is the first derivative matrix.
//' @param v_hat A \code{mat} that is the bootstrapped V matrix
//' @param omega A \code{mat} that is the inverse of the diagonal V matrix.
//' @param alpha A \code{double} that contains the confidence level.
//' @return A \code{mat} that has the first column 
//' @backref src/inference.cpp
//' @backref src/inference.h
//' @keywords internal
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
//' @template tsobj_cpp
//' @param model_type A \code{string} that contains the model type: \code{"imu"} or \code{"ssm"}
//' @param tau A \code{vec} containing the scales of a proccess.
//' @param v_hat A \code{mat} that contains the bootstrapped matrix.
//' @param wv_empir A \code{vec} that contains the empirical wavelet variance.
//' @return A \code{vec} that has
//' \itemize{
//' \item Test Statistic
//' \item P-Value
//' \item DF
//' } 
//' @backref src/inference.cpp
//' @backref src/inference.h
//' @keywords internal
// [[Rcpp::export]]
arma::vec gof_test(arma::vec theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   std::string model_type,
                   const arma::vec& tau,
                   const arma::mat& v_hat, const arma::vec& wv_empir){
  
  arma::mat omega = arma::inv(v_hat);
  
  arma::vec estimate = gmwm_engine(theta,
                                   desc, objdesc, 
                                  model_type, 
                                  wv_empir,
                                  v_hat,
                                  tau,
                                  false); // starting is false
  
  

  double test_stat = getObjFun(estimate, desc, objdesc, model_type,
                               omega, wv_empir, tau);
    
  unsigned int df = tau.n_elem - theta.n_elem;
  
  double p_value = 1.0 - R::pchisq(test_stat, df, true, false);

  arma::vec out(3);
  out(0) = test_stat;
  out(1) = p_value;
  out(2) = df;

  return out;
}

//' @title Compute the Bootstrapped GoF Test
//' @description Handles the bootstrap computation and the bootstrapped p-value.
//' @param obj_value A \code{double} that contains the optimized objective function value.
//' @param bs_obj_values A \code{vec} that contains the objective function values under bootstrap.
//' @param alpha A \code{double} that indicates the confidence.
//' @param bs_gof_p_ci A \code{bool} that indicates whether CIs should be included or not.
//' @return A \code{vec} that has
//' \itemize{
//' \item Test Statistic
//' \item Low CI
//' \item Upper CI - BS
//' } 
//' @backref src/inference.cpp
//' @backref src/inference.h
//' @keywords internal
// [[Rcpp::export]]
arma::vec bootstrap_gof_test(double obj_value, arma::vec bs_obj_values, double alpha, bool bs_gof_p_ci){
  
  arma::vec temp(1+2*bs_gof_p_ci);
  temp(0) = sum(obj_value < bs_obj_values)/double(bs_obj_values.n_elem);
  if(bs_gof_p_ci){
    temp.rows(1,2) = boot_pval_gof(obj_value, bs_obj_values, 1000, alpha);
  }
  return temp;
}
