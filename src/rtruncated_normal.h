#ifndef RTRUNCATED_NORMAL_H
#define RTRUNCATED_NORMAL_H

// One
arma::vec rtruncated_normal(unsigned int n, double mu, double sigma, double a, double b);
double rtruncated_normal(double mu, double sigma, double a, double b);

double sim_truncated_normal(double phi_a_cdf, double phi_b_cdf, double mu, double sigma);

#endif