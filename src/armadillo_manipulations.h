#ifndef ARMADILLO_MANIPULATIONS
#define ARMADILLO_MANIPULATIONS

arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end);

arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end);

arma::vec reverse_vec(arma::vec x);

arma::mat field_to_matrix(arma::field<arma::vec> x);

double sum_field_vec(const arma::field<arma::vec>& x);

#endif
