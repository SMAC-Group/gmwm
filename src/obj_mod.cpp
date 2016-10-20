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


#include <RcppArmadillo.h>


//' @title Extract Object
//' @description Extracts the object information and returns it.
//' @param theta A \code{vec} containing the theta values.
//' @param objdesc A \code{vec} at the desc point.
//' @param cur_position An \code{integer} at the current position.
//' @return A \code{field<vec>} containing the breakdown of the object.
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::vec> obj_extract(arma::vec theta,
                                   arma::vec objdesc,
                                   unsigned int& cur_position) {
  
  
  unsigned int nobs = objdesc.n_elem;
  
  arma::field<arma::vec> out(nobs);
  
  unsigned int i;
  
  for(i = 0; i < nobs; i++){
    unsigned int obj = objdesc(i);
    out(i) = theta.rows(cur_position, cur_position + obj - 1);
    cur_position += obj;
  }
  
  return out;
}