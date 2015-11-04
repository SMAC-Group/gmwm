/* Copyright (C) 2014 - 2015  James Balamuta
 *
 * This file is part of GMWM R Methods Package
 *
 * The file uses methods in the r-to-armadillo project and is free software: you can redistribute it and/or modify it
 * under the terms of the MIT License.
 *
 * The r-to-armadillo project is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

#ifndef ARMADILLO_MANIPULATIONS
#define ARMADILLO_MANIPULATIONS

arma::mat sort_mat(arma::mat x, unsigned int col);
  
arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end);

arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end);

arma::vec reverse_vec(arma::vec x);

arma::mat field_to_matrix(arma::field<arma::vec> x);

double sum_field_vec(const arma::field<arma::vec>& x);

#endif
