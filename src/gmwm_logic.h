#ifndef GMWM_FUNCTIONS
#define GMWM_FUNCTIONS

arma::vec gmwm_cpp(const arma::vec& theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type, 
                          const arma::mat& omega, const arma::vec& wv_empir,
                          const arma::vec& tau);
                          
arma::vec adv_gmwm_cpp(const arma::vec& theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type, 
                          const arma::mat& omega, const arma::vec& wv_empir,
                          const arma::vec& tau);
                          
arma::mat gmwm_bootstrapper(const arma::vec&  theta,
                            const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                            unsigned int N, bool robust, double eff,
                            unsigned int H = 100, bool diagonal_matrix = true);
                            
                          
arma::vec gmwm_engine(const arma::vec& theta,
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                      std::string model_type,
                      arma::vec wv_empir,
                      arma::mat omega,
                      arma::vec scales,
                      bool starting = true);

arma::field<arma::mat> gmwm_update_cpp(arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, unsigned int N, double expect_diff, 
                                      arma::mat orgV, arma::vec scales, arma::vec wv_empir,
                                      bool starting = true, 
                                      std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                      unsigned int G = 1000, 
                                      bool robust=false, double eff = 0.6, bool inference = false);
                                      
arma::field<arma::mat> gmwm_master_cpp(const arma::vec& data, 
                                      arma::vec theta,
                                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, 
                                      std::string model_type, bool starting = true,
                                      double alpha = 0.05, 
                                      std::string compute_v = "fast", unsigned int K = 1, unsigned int H = 100,
                                      unsigned int G = 1000, 
                                      bool robust=false, double eff = 0.6, bool inference = false, bool modelselect = false);
                                      
#endif
