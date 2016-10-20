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

#ifndef INFERENCE_H
#define INFERENCE_H

arma::mat calculate_psi_matrix(const arma::mat& A, const arma::mat& v_hat, const arma::mat& omega);

arma::mat format_ci(const arma::vec& theta,
                    const arma::vec& se,
                    double z);

arma::mat theta_ci(const arma::vec& theta,
                   const arma::mat& A, 
                   const arma::mat& v_hat, const arma::mat& omega, double alpha);

arma::vec gof_test(arma::vec theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   std::string model_type,
                   const arma::vec& tau,
                   const arma::mat& v_hat, const arma::vec& wv_empir);
  
arma::vec bootstrap_gof_test(double obj_value, arma::vec bs_obj_values, double alpha, bool bs_gof_p_ci);
#endif
