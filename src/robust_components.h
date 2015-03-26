#ifndef ROBUST_COMPONENTS
#define ROBUST_COMPONENTS

double objFun_find_biwc(double crob, double eff);

double find_biwc(double eff);

double objFun_sig_rob_bw(double sig2_bw, arma::vec x, double a_of_c, double crob_bw);

double sig_rob_bw(arma::vec y, double eff);

#endif