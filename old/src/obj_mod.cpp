#include <RcppArmadillo.h>


//' @title Extract Object
//' @description Extracts the object information and returns it.
//' @param theta A \code{vec} containing the theta values.
//' @param objdesc A \code{vec} at the desc point.
//' @param cur_position An \code{integer} at the current position.
//' @return A \code{field<vec>} containing the breakdown of the object.
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::vec> obj_extract(arma::vec theta,
                                   arma::vec objdesc,
                                   unsigned int& cur_position) {
  
  
  unsigned int nobs = objdesc.n_elem;
  
  arma::field<arma::vec> out(nobs);
  
  unsigned int i;
  
  for(i = 0; i < nobs; i++){
    unsigned int obj = objdesc(i);
    out(i) = theta.rows(cur_position, cur_position + obj - 1);
    cur_position += obj;
  }
  
  return out;
}