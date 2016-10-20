/* Copyright (C) 2016  James Balamuta
 *
 * This file is part of GMWM R Methods Package
 *
 * The file uses methods in the rgen project and is free software: you can redistribute it and/or modify it
 * under the terms of the MIT License.
 *
 * The rgen project is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

#ifndef RTRUNCATED_NORMAL_H
#define RTRUNCATED_NORMAL_H

// Multiple obs
arma::vec rtruncated_normal(unsigned int n, double mu, double sigma, double a, double b);

// One obs
double rtruncated_normal(double mu, double sigma, double a, double b);

// RNG Sampler
double sim_truncated_normal(double phi_a_cdf, double phi_b_cdf, double mu, double sigma);

#endif