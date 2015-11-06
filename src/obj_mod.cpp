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


#include <RcppArmadillo.h>
using namespace Rcpp;


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