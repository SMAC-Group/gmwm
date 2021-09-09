#ifndef GMWM_FUNCTIONS
#define GMWM_FUNCTIONS

arma::vec gmwm_engine(const arma::vec& theta,
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                      std::string model_type,
                      arma::vec wv_empir,
                      arma::mat omega,
                      arma::vec scales,
                      bool starting = true);

arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                       const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                       std::string model_type, unsigned int N, double expect_diff, double ranged, 
                                       const arma::mat& orgV, const arma::vec& scales, const arma::mat& wv,
                                       bool starting = true, 
                                       std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                       unsigned int G = 1000, 
                                       bool robust=false, double eff = 0.6);
                                      
arma::field<arma::mat> gmwm_master_cpp(arma::vec& data, 
                                       arma::vec theta,
                                       const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                       std::string model_type, bool starting = true,
                                       double alpha = 0.05, 
                                       std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                       unsigned int G = 1000, 
                                       bool robust=false, double eff = 0.6);

arma::field<arma::mat> gmwm_master_wv_cpp(arma::mat wvar,
                                          unsigned int N,
                                          double expect_diff,
                                          arma::mat omega,
                                          double ranged,
                                          arma::vec theta,
                                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                          std::string model_type, bool starting = true,
                                          double alpha = 0.05, 
                                          std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                          unsigned int G = 1000, 
                                          bool robust=false, double eff = 0.6);
  
                                      
#endif
