#ifndef INFERENCE_H
#define INFERENCE_H

arma::mat calculate_psi_matrix(const arma::mat& A, const arma::mat& v_hat, const arma::mat& omega);

arma::mat format_ci(const arma::vec& theta,
                    const arma::vec& se,
                    double z);

arma::mat theta_ci(const arma::vec& theta,
                   const arma::mat& A, 
                   const arma::mat& v_hat, const arma::mat& omega, double alpha);

arma::vec gof_test(arma::vec theta, 
                   const std::vector<std::string>& desc,
                   const arma::field<arma::vec>& objdesc,
                   std::string model_type,
                   const arma::vec& tau,
                   const arma::mat& v_hat, const arma::vec& wv_empir);
  
arma::vec bootstrap_gof_test(double obj_value, arma::vec bs_obj_values, double alpha, bool bs_gof_p_ci);
#endif
