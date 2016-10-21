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

#ifndef GEN_PROCESS
#define GEN_PROCESS

arma::vec gen_wn(const unsigned int N, const double sigma2);

arma::vec gen_dr(const unsigned int N, const double omega);

arma::vec gen_qn(const unsigned int N, double q2);

arma::vec gen_ar1(const unsigned int N, const double phi, const double sigma2);

arma::vec gen_rw(const unsigned int N, const double sigma2);

arma::vec gen_arma(const unsigned int N,
                   const arma::vec& ar, const arma::vec& ma,
                   const double sigma2, 
                   unsigned int n_start);

arma::vec gen_arima(const unsigned int N,
                    const arma::vec& ar,
                    const unsigned int d,
                    const arma::vec& ma,
                    const double sigma2, 
                    unsigned int n_start);

arma::vec gen_sarma(const unsigned int N,
                    const arma::vec& ar, const arma::vec& ma,
                    const arma::vec& sar, const arma::vec& sma,
                    const double sigma2, 
                    unsigned int s, 
                    unsigned int n_start);

arma::vec gen_generic_sarima(const unsigned int N,
                             const arma::vec& theta_values, 
                             const arma::vec& objdesc,
                             double sigma2,
                             unsigned int n_start);

arma::vec gen_model(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc);

arma::mat gen_lts_cpp(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc);

#endif
