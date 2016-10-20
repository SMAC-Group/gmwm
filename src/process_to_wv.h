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

#ifndef PROCESS_TO_WV
#define PROCESS_TO_WV

arma::vec arma_to_wv(arma::vec ar, arma::vec ma, double sigma2, arma::vec tau);

arma::vec arma11_to_wv(double phi, double theta, double sigma2, const arma::vec& tau);

arma::vec ar1_to_wv(double phi, double sigma2, const arma::vec& tau);

arma::vec ma1_to_wv(double theta, double sigma2, const arma::vec& tau);

arma::vec qn_to_wv(double q2, const arma::vec& tau);

arma::vec wn_to_wv(double sigma2, arma::vec tau);

arma::vec rw_to_wv(double gamma2, const arma::vec& tau);

arma::vec dr_to_wv(double omega, const arma::vec& tau);

arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau);
                         
arma::mat decomp_theoretical_wv(const arma::vec& theta, 
                                const std::vector<std::string>& desc,
                                const arma::field<arma::vec>& objdesc, const arma::vec& tau);

arma::vec decomp_to_theo_wv(const arma::mat& decomp);

#endif
