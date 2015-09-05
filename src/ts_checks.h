#ifndef TS_CHECKS_H
#define TS_CHECKS_H

double minroot(const arma::cx_vec& x);

bool invert_check(const arma::vec& x);

std::map<std::string, int> count_models(const std::vector<std::string>& desc);

arma::vec order_AR1s(arma::vec theta, const std::vector<std::string>& desc, const arma::field<arma::vec> objdesc);

#endif