#ifndef PROCESS_TO_WV
#define PROCESS_TO_WV

arma::vec arma_to_wv(arma::vec ar, arma::vec ma, arma::vec tau, double sigma);

arma::vec qn_to_wv(double q2, const arma::vec& Tau);

arma::vec wn_to_wv(double sig2, arma::vec Tau);

arma::vec rw_to_wv(double sig2, const arma::vec& Tau);

arma::vec dr_to_wv(double omega,const arma::vec& Tau);

arma::vec ar1_to_wv(double phi, double sig2, const arma::vec& Tau);

arma::vec theoretical_wv(const arma::vec& theta, 
                         const std::vector<std::string>& desc,
                         const arma::field<arma::vec>& objdesc, const arma::vec& tau);
#endif
