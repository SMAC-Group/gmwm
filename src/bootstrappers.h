#ifndef BOOTSTRAPPER_H
#define BOOTSTRAPPER_H

arma::mat cov_bootstrapper(const arma::vec&  theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                           unsigned int N, bool robust, double eff,
                           unsigned int H, bool diagonal_matrix);

arma::vec gmwm_sd_bootstrapper(const arma::vec&  theta,
                               const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                               const arma::vec& scales, std::string model_type,
                               unsigned int N, bool robust, double eff, double alpha,
                               unsigned int H);

arma::mat optimism_bootstrapper(const arma::vec&  theta,
                                const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                const arma::vec& scales, std::string model_type, 
                                unsigned int N, bool robust, double eff, double alpha,
                                unsigned int H);

arma::vec boot_pval_gof(double obj, const arma::vec& obj_boot, unsigned int B, double alpha);

arma::field<arma::mat> all_bootstrapper(const arma::vec&  theta,
                                        const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                                        const arma::vec& scales, std::string model_type, 
                                        unsigned int N, bool robust, double eff, double alpha,
                                        unsigned int H);

#endif