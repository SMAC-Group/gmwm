#ifndef MODEL_SELECTION
#define MODEL_SELECTION

arma::mat D_matrix(arma::mat At_j, arma::mat omega, arma::mat diff);

arma::mat B_matrix(arma::mat A, arma::mat a_omega);

double model_score(arma::mat A, arma::mat At_j, arma::mat omega, arma::mat v_hat, arma::vec diff, unsigned int T);

#endif
