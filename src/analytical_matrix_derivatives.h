/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
 * (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#ifndef ANALYTICAL_MATRIX_DERIVATIVES
#define ANALYTICAL_MATRIX_DERIVATIVES

arma::mat deriv_AR1(double phi, double sig2, arma::vec tau);

arma::mat deriv_2nd_AR1(double phi, double sig2, arma::vec tau);

arma::mat deriv_DR(double omega, arma::vec tau);

arma::mat deriv_2nd_DR(arma::vec tau);

arma::mat deriv_QN(arma::vec tau);

arma::mat deriv_RW(arma::vec tau);

arma::mat deriv_WN(arma::vec tau);

arma::mat derivative_first_matrix(const arma::vec& theta, 
                                  const std::vector<std::string>& desc,
                                  const arma::field<arma::vec>& objdesc,
                                  const arma::vec& tau);

arma::mat D_matrix(const arma::vec& theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   const arma::vec& tau, const arma::vec& omegadiff);

#endif
