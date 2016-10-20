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

#ifndef BOOTSTRAPPER_H
#define BOOTSTRAPPER_H

arma::mat cov_bootstrapper(const arma::vec&  theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                           unsigned int N, bool robust, double eff,
                           unsigned int H, bool diagonal_matrix);

arma::vec gmwm_sd_bootstrapper(const arma::vec&  theta,
                               const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                               const arma::vec& scales, std::string model_type,
                               unsigned int N, bool robust, double eff, double alpha,
                               unsigned int H);

arma::mat optimism_bootstrapper(const arma::vec&  theta,
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                const arma::vec& scales, std::string model_type, 
                                unsigned int N, bool robust, double eff, double alpha,
                                unsigned int H);

arma::field<arma::mat> opt_n_gof_bootstrapper(const arma::vec&  theta,
                                              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                              const arma::vec& scales, std::string model_type, 
                                              unsigned int N, bool robust, double eff, double alpha,
                                              unsigned int H);

arma::vec boot_pval_gof(double obj, const arma::vec& obj_boot, unsigned int B, double alpha);

arma::field<arma::mat> all_bootstrapper(const arma::vec&  theta,
                                        const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                        const arma::vec& scales, std::string model_type, 
                                        unsigned int N, bool robust, double eff, double alpha,
                                        unsigned int H);

#endif