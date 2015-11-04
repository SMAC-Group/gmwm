/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Q Public License included within the packages source
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Q Public License
 * along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#ifndef POLYROOT_H
#define POLYROOT_H

double myfmod_cpp(double x1, double x2);

double R_pow_cpp(double x, double y);

double R_pow_di_cpp(double x, int n);


// Function definitions
void calct_cpp(bool &bol);
bool fxshft_cpp(const int l2, double &zr, double &zi);
bool vrshft_cpp(const int l3, double &zr, double &zi);
void nexth_cpp(const bool bol);
void noshft_cpp(const int l1);


void polyev_cpp(const int n, const double s_r, const double s_i, 
                       std::vector<double> &p_r, std::vector<double> &p_i, 
                       std::vector<double> &q_r, std::vector<double> &q_i, 
                       double &v_r, double &v_i);
double errev_cpp(const int n, 
                        const std::vector<double> &qr, std::vector<double> &qi,
                        const double ms, const double mp, const double a_re, const double m_re);
double cpoly_cauchy_cpp(const int n, std::vector<double> &pot, std::vector<double> &q);
double cpoly_scale_cpp(const int n, std::vector<double> &pot,
                              const double eps, const double BIG,
                              const double small, const double base);
void cdivid_cpp(const double ar, const double ai, 
                       const double br, const double bi, 
                       double &cr, double &ci);

void polyroot_cpp(const std::vector<double> &opr, const std::vector<double> &opi, int &degree,
                         std::vector<double> &zeror, std::vector<double> &zeroi, bool &fail);

arma::cx_vec do_polyroot_arma(const arma::cx_vec& z);

std::vector< std::complex<double> > do_polyroot_cpp(const std::vector< std::complex<double> >& z);

#endif