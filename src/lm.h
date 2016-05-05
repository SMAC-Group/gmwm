#ifndef LM_H
#define LM_H

arma::field<arma::vec> lm_arma(const arma::vec & y, const arma::mat & X);
arma::field<arma::vec> lm_dr(const arma::vec & x);

#endif