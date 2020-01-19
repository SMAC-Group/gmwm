#ifndef SARMA_H
#define SARMA_H

arma::vec sarma_objdesc(const arma::vec& ar, const arma::vec&ma,
                        const arma::vec& sar, const arma::vec& sma,
                        int s = 12, int i = 0, int si = 0);

arma::vec sarma_components(const arma::vec& objdesc);


arma::field<arma::vec> sarma_expand(const arma::vec& params, const arma::vec& objdesc);

arma::field<arma::vec> sarma_expand_unguided(const arma::vec& params,
                                             unsigned int np, unsigned int nq,
                                             unsigned int nsp, unsigned int nsq,
                                             unsigned int ns,
                                             unsigned int p,
                                             unsigned int q);

arma::vec sarma_params_construct(const arma::vec& ar, const arma::vec& ma,
                                 const arma::vec& sar, const arma::vec& sma);

arma::vec sarma_calculate_spadding(unsigned int np, unsigned int nq,
                                   unsigned int nsp, unsigned int nsq,
                                   unsigned int ns);
#endif