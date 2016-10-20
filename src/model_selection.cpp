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

#include "model_selection.h"

// Needed for COV bootstrap
#include "bootstrappers.h"


/* To whom it may concern....
 * 
 * The D matrix is located in the analytical_matrix_derivatives.cpp code file.
 * 
 */

//' @title B Matrix
//' @description B Matrix
//' @param A A \code{mat} containing the first derivatives of the process.
//' @param at_omega A \code{mat} containing A^T * Omega
//' @return A \code{mat}
//' @backref src/model_selection.cpp
//' @backref src/model_selection.h
//' @keywords internal
// [[Rcpp::export]]
arma::mat B_matrix(const arma::mat& A, const arma::mat& at_omega){
  unsigned int p = A.n_cols;
  arma::mat B(p,p);
  
  for(unsigned int i = 0; i < p; i++){
    B.col(i) = at_omega* A.col(i);
  }
  
  return B;
}

//' @title Model Score
//' @description Calculates the modeling score of a GMWM
//' @param A A \code{mat} that contains the first derivatives of the processes
//' @param At_j  A \code{mat} that contains the second derivative of each process
//' @param omega A \code{mat} that contains the omega used when calculating the GMWM
//' @param v_hat A \code{mat} that contains the covariance matrix
//' @param diff A \code{vec} that is the difference of the WV empirical and WV theoretical
//' @return A \code{vec}
//' @keywords internal
//' @backref src/model_selection.cpp
//' @backref src/model_selection.h
//' @details
//' The equation is slightly different than that stated in the paper due to the bootstrap already incorporating in 
//' N.
// [[Rcpp::export]]
arma::vec model_score(arma::mat A, arma::mat D, arma::mat omega, arma::mat v_hat, double obj_value){
  
  arma::mat At = arma::trans(A);
  arma::mat B = B_matrix(A, At*omega);
  
  arma::mat d_b = D-B;

  arma::mat db_t = arma::trans(d_b);
  
  arma::mat dTheta = -1*arma::pinv(db_t * d_b)*db_t*At*omega;

  arma::vec score_info(2);
  
  double optimism = 2.0*arma::as_scalar(arma::trace(A * dTheta * omega * v_hat));
  
  score_info(0) = obj_value + optimism;
  score_info(1) = optimism;
  
  // Return the score
  return score_info;
}

// arma::vec mscore(){
//   
//   
//     arma::vec score;
//     /* A note to someone in the future...
//     * Yes, there is a difference in order between the diff (wv_empir-theo) for D_matrix
//     *  and the model_score diff (theo-wv_empir).
//     */
//     if(bs_optimism){
//       arma::vec temp(2);
//       
//       double optimism = 2*sum(diagvec(cov_nu_nu_theta * omega));
//       
//       temp(0) = obj_value + optimism;
//       temp(1) = optimism;
//       
//       score = temp;
//     }else{
//       // Create the D Matrix (note this is in the analytical_matrix_derivaties.cpp file)
//       arma::mat D = D_matrix(theta, desc, objdesc, scales, omega*(wv_empir - theo));
//       
//       // Calculate the model score according to model selection criteria paper
//       score = model_score(A, D, omega, V,  obj_value);
//     }
// 
// }
