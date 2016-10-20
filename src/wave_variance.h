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

#ifndef WAVE_VARIANCE
#define WAVE_VARIANCE

arma::mat ci_eta3(const arma::vec& y, const arma::vec& dims, double alpha_ov_2);

arma::mat ci_eta3_robust(const arma::vec& wv_robust, const arma::mat& wv_ci_class, double alpha_ov_2, double eff);

arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& wv, 
                            std::string type, double alpha_ov_2, bool robust, double eff);
 
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw,
                        bool robust, double eff);

arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt_bw,
                     bool robust=false, double eff=0.6, double alpha = 0.05, 
                     std::string ci_type="eta3");

arma::mat modwt_wvar_cpp(const arma::vec& signal, unsigned int nlevels = 4,
                         bool robust=false, double eff=0.6, double alpha = 0.05, 
                         std::string ci_type="eta3", std::string strWavelet="haar", std::string decomp="modwt");

arma::field<arma::mat> batch_modwt_wvar_cpp(const arma::mat& signal, unsigned int nlevels = 4,
                                            bool robust=false, double eff=0.6, double alpha = 0.05, 
                                            std::string ci_type="eta3", std::string strWavelet="haar", std::string decomp="modwt");

arma::vec scales_cpp(unsigned int nb_level);

#endif
