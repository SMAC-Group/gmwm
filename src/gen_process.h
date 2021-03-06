#ifndef GEN_PROCESS
#define GEN_PROCESS

arma::vec gen_wn(const unsigned int N, const double sigma2);

arma::vec gen_dr(const unsigned int N, const double omega);

arma::vec gen_qn(const unsigned int N, double q2);

arma::vec gen_ar1(const unsigned int N, const double phi, const double sigma2);

arma::vec gen_rw(const unsigned int N, const double sigma2);

arma::vec gen_arma(const unsigned int N,
                   const arma::vec& ar, const arma::vec& ma,
                   const double sigma2, 
                   unsigned int n_start);

arma::vec gen_arima(const unsigned int N,
                    const arma::vec& ar,
                    const unsigned int d,
                    const arma::vec& ma,
                    const double sigma2, 
                    unsigned int n_start);

arma::vec gen_sarma(const unsigned int N,
                    const arma::vec& ar, const arma::vec& ma,
                    const arma::vec& sar, const arma::vec& sma,
                    const double sigma2, 
                    unsigned int s, 
                    unsigned int n_start);

arma::vec gen_generic_sarima(const unsigned int N,
                             const arma::vec& theta_values, 
                             const arma::vec& objdesc,
                             double sigma2,
                             unsigned int n_start);

arma::vec gen_model(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc);

arma::mat gen_lts_cpp(unsigned int N, const arma::vec& theta, const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc);

#endif
