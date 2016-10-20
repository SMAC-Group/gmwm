/* Copyright (C) 2016 James Balamuta, Stephane Guerrier, Roberto Molinari
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

#ifndef SARMA_H
#define SARMA_H

arma::vec sarma_objdesc(const arma::vec& ar, const arma::vec&ma,
                        const arma::vec& sar, const arma::vec& sma,
                        int s = 12, int i = 0, int si = 0);

arma::vec sarma_components(const arma::vec& objdesc);


arma::field<arma::vec> sarma_expand(const arma::vec& params, const arma::vec& objdesc);

arma::field<arma::vec> sarma_expand_unguided(const arma::vec& params,
                                             unsigned int np, unsigned int nq,
                                             unsigned int nsp, unsigned int nsq,
                                             unsigned int ns,
                                             unsigned int p,
                                             unsigned int q);

arma::vec sarma_params_construct(const arma::vec& ar, const arma::vec& ma,
                                 const arma::vec& sar, const arma::vec& sma);

arma::vec sarma_calculate_spadding(unsigned int np, unsigned int nq,
                                   unsigned int nsp, unsigned int nsq,
                                   unsigned int ns);
#endif