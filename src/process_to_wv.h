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

#ifndef PROCESS_TO_WV
#define PROCESS_TO_WV

arma::vec arma_to_wv(arma::vec ar, arma::vec ma, arma::vec tau, double sigma);

arma::vec qn_to_wv(double q2, const arma::vec& Tau);

arma::vec wn_to_wv(double sig2, arma::vec Tau);

arma::vec rw_to_wv(double sig2, const arma::vec& Tau);

arma::vec dr_to_wv(double omega,const arma::vec& Tau);

arma::vec ar1_to_wv(double phi, double sig2, const arma::vec& Tau);

arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau);
                         
arma::mat decomp_theoretical_wv(const arma::vec& theta, 
                                const std::vector<std::string>& desc,
                                const arma::field<arma::vec>& objdesc, const arma::vec& tau);

arma::vec decomp_to_theo_wv(const arma::mat& decomp);

#endif
