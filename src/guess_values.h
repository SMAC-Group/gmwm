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

#ifndef GUESS_VALUES
#define GUESS_VALUES
#include <map>

arma::vec ar1_draw(unsigned int draw_id, double last_phi, double sigma_tot, std::string model_type);

arma::vec arma_draws(unsigned int p, unsigned int q, double sigma2_total);

double dr_slope(const arma::vec& data);

std::string dom_process(double first_wv, double ci_low, double ci_high);

arma::vec guess_initial(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::mat& wv, const arma::vec& tau, double ranged, unsigned int G=1000);

arma::vec guess_initial_old(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                            std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                            const arma::vec& wv_empir, const arma::vec& tau, unsigned int B = 1000);


#endif
