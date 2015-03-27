#ifndef COVARIANCE_MATRIX
#define COVARIANCE_MATRIX

arma::field<arma::mat> compute_cov_cpp(arma::field<arma::vec> signal_modwt, unsigned int nb_level, std::string compute_v,
                                        bool robust, double eff);

#endif
