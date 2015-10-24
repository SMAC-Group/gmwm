#ifndef WAVE_VARIANCE
#define WAVE_VARIANCE

arma::mat ci_eta3(arma::vec y,  arma::vec dims, double alpha_ov_2);

arma::mat ci_eta3_robust(arma::vec wv_robust, arma::mat wv_ci_class, double alpha_ov_2, double eff);

arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& wv, 
                            std::string type, double alpha_ov_2, bool robust, double eff);
 
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw,
                        bool robust, double eff);

arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt, 
                   bool robust, double eff, double alpha, 
                   std::string ci_type, std::string strWavelet);

arma::vec scales_cpp(unsigned int nb_level);

#endif
