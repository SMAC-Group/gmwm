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

#ifndef TRANSFORM_DATA
#define TRANSFORM_DATA

arma::vec pseudo_logit_inv(const arma::vec& x);

arma::vec logit_inv(const arma::vec& x);

arma::vec pseudo_logit(const arma::vec& x);

double pseudo_logit(double x);

arma::vec logit(const arma::vec& x);

double logit(double x);

arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type);
                           
arma::vec untransform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type);                           

#endif
