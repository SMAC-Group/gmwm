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

#ifndef R_TO_ARMADILLO
#define R_TO_ARMADILLO

#define my_isok(x) (!ISNA(x) & !ISNAN(x))

arma::vec seq_cpp(int a, int b);

arma::vec seq_len_cpp(unsigned int n);

arma::vec quantile_cpp(arma::vec x, const arma::vec& probs);

arma::vec diff_cpp(arma::vec x, unsigned int lag = 1, unsigned int differences = 1);

arma::vec ARMAtoMA_cpp(arma::vec ar, arma::vec ma, int lag_max);

arma::vec cfilter(arma::vec x, arma::vec filter, int sides = 2, bool circular = false);

arma::vec rfilter(arma::vec x, arma::vec filter, arma::vec init);

arma::mat expand_grid_red(int nx);

arma::vec ARMAacf_cpp(arma::vec ar, arma::vec ma, unsigned int lag_max);

arma::vec dft_acf(const arma::vec& x);

double mean_diff(const arma::vec& x);

arma::vec num_rep(const arma::vec& x, unsigned int n);

arma::vec intgr_vec(const arma::vec& x, const arma::vec& xi, unsigned int lag);

arma::vec diff_inv_values(const arma::vec& x, unsigned int lag, unsigned int d, const arma::vec& xi);

arma::vec diff_inv(const arma::vec& x, unsigned int lag, unsigned int d);

#endif
