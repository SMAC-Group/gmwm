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

#ifndef GMWM_FUNCTIONS
#define GMWM_FUNCTIONS

arma::mat cov_bootstrapper(const arma::vec&  theta,
                            const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                            unsigned int N, bool robust, double eff,
                            unsigned int H = 100, bool diagonal_matrix = true);
 
 arma::vec gmwm_sd_bootstrapper(const arma::vec&  theta,
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                const arma::vec& scales, std::string model_type,
                                unsigned int N, bool robust, double eff, double alpha,
                                unsigned int H);                           
                          
arma::vec gmwm_engine(const arma::vec& theta,
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                      std::string model_type,
                      arma::vec wv_empir,
                      arma::mat omega,
                      arma::vec scales,
                      bool starting = true);

arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, unsigned int N, double expect_diff, 
                                      arma::mat orgV, arma::vec scales, arma::vec wv_empir,
                                      bool starting = true, 
                                      std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                      unsigned int G = 1000, 
                                      bool robust=false, double eff = 0.6);
                                      
arma::field<arma::mat> gmwm_master_cpp(const arma::vec& data, 
                                      arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, bool starting = true,
                                      double alpha = 0.05, 
                                      std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                      unsigned int G = 1000, 
                                      bool robust=false, double eff = 0.6);
                                      
#endif
