#ifndef ANALYTICAL_MATRIX_DERIVATIVES
#define ANALYTICAL_MATRIX_DERIVATIVES

arma::mat deriv_AR1(double phi, double sig2, arma::vec tau);

arma::mat deriv_2nd_AR1(double phi, double sig2, arma::vec tau);

arma::mat deriv_DR(double omega, arma::vec tau);

arma::mat deriv_2nd_DR(arma::vec tau);

arma::mat deriv_QN(arma::vec tau);

arma::mat deriv_RW(arma::vec tau);

arma::mat deriv_WN(arma::vec tau);

#endif
