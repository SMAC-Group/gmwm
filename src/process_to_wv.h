#ifndef PROCESS_TO_WV
#define PROCESS_TO_WV

arma::vec arma_to_wv(arma::vec ar, arma::vec ma, double sigma2, arma::vec tau);

arma::vec arma11_to_wv(double phi, double theta, double sigma2, const arma::vec& tau);

arma::vec ar1_to_wv(double phi, double sigma2, const arma::vec& tau);

arma::vec ma1_to_wv(double theta, double sigma2, const arma::vec& tau);

arma::vec qn_to_wv(double q2, const arma::vec& tau);

arma::vec wn_to_wv(double sigma2, arma::vec tau);

arma::vec rw_to_wv(double gamma2, const arma::vec& tau);

arma::vec dr_to_wv(double omega, const arma::vec& tau);

arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau);
                         
arma::mat decomp_theoretical_wv(const arma::vec& theta, 
                                const std::vector<std::string>& desc,
                                const arma::field<arma::vec>& objdesc, const arma::vec& tau);

arma::vec decomp_to_theo_wv(const arma::mat& decomp);

#endif
