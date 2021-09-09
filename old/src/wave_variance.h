#ifndef WAVE_VARIANCE
#define WAVE_VARIANCE

arma::mat ci_eta3(const arma::vec& y, const arma::vec& dims, double alpha_ov_2);

arma::mat ci_eta3_robust(const arma::vec& wv_robust, const arma::mat& wv_ci_class, double alpha_ov_2, double eff);

arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& wv, 
                            std::string type, double alpha_ov_2, bool robust, double eff);
 
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw,
                        bool robust, double eff);

arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt_bw,
                     bool robust=false, double eff=0.6, double alpha = 0.05, 
                     std::string ci_type="eta3");

arma::mat modwt_wvar_cpp(const arma::vec& signal, unsigned int nlevels = 4,
                         bool robust=false, double eff=0.6, double alpha = 0.05, 
                         std::string ci_type="eta3", std::string strWavelet="haar", std::string decomp="modwt");

arma::field<arma::mat> batch_modwt_wvar_cpp(const arma::mat& signal, unsigned int nlevels = 4,
                                            bool robust=false, double eff=0.6, double alpha = 0.05, 
                                            std::string ci_type="eta3", std::string strWavelet="haar", std::string decomp="modwt");

arma::vec scales_cpp(unsigned int nb_level);

#endif
