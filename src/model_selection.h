#ifndef MODEL_SELECTION
#define MODEL_SELECTION

arma::mat B_matrix(const arma::mat& A, const arma::mat& at_omega);

arma::vec model_score(arma::mat A, arma::mat D, arma::mat omega, arma::mat v_hat, arma::vec diff);

#endif
