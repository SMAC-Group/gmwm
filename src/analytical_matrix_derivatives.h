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

#ifndef ANALYTICAL_MATRIX_DERIVATIVES
#define ANALYTICAL_MATRIX_DERIVATIVES

arma::mat deriv_AR1(double phi, double sigma2, const arma::vec& tau);

arma::mat deriv_2nd_AR1(double phi, double sigma2, const arma::vec& tau);

arma::mat deriv_MA1(double theta, double sigma2, const arma::vec& tau);

arma::mat deriv_2nd_MA1(double theta, double sigma2, const arma::vec& tau);

arma::mat deriv_DR(double omega, const arma::vec& tau);

arma::mat deriv_2nd_DR(const arma::vec& tau);

arma::mat deriv_QN(const arma::vec& tau);

arma::mat deriv_RW(const arma::vec& tau);

arma::mat deriv_WN(const arma::vec& tau);

arma::mat derivative_first_matrix(const arma::vec& theta, 
                                  const std::vector<std::string>& desc,
                                  const arma::field<arma::vec>& objdesc,
                                  const arma::vec& tau);

arma::mat D_matrix(const arma::vec& theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   const arma::vec& tau, const arma::vec& omegadiff);

#endif
