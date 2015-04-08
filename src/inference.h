#ifndef INFERENCE
#define INFERENCE

arma::mat calculate_psi_matrix(const arma::mat& D, const arma::mat& v_hat, const arma::mat& omega);


arma::mat theta_ci(arma::vec theta, arma::mat psi, double alpha = 0.05);


arma::vec gof_test(const arma::vec& theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   std::string model_type,
                   const arma::vec& tau,
                   const arma::mat v_hat, 
                   const arma::vec& wv_empir);
  

arma::field<arma::mat> inference_summary(const arma::vec& theta, 
                                        const std::vector<std::string>& desc,
                                        const arma::field<arma::vec>& objdesc,
                                        std::string model_type,
                                        const arma::vec& tau,
                                        arma::mat v_hat, arma::mat omega, arma::vec wv_empir, double alpha = 0.05);

#endif
