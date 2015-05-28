#ifndef TS_CHECKS_H
#define TS_CHECKS_H

double minroot(const arma::cx_vec& x);

bool invert_check(const arma::cx_vec& x);

std::map<std::string, int> count_models(const std::vector<std::string>& desc);

#endif