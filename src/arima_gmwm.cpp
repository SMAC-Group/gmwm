#include <RcppArmadillo.h>

#include "arima_gmwm.h"

using namespace Rcpp;


// [[Rcpp::export]]
arma::vec Rcpp_ARIMA(const arma::vec& data,
                     const arma::vec& params){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function arima = stats["arima"];    
  
  arma::vec aparams(3);
  
  aparams(0) = params(0);
  aparams(1) = 0;
  aparams(2) = params(1);
  
  Rcpp::List Opt= arima(_["x"] = data,
                        _["order"] = aparams,
                        _["include.mean"] = false,
                        _["method"] = "CSS");
  
  arma::vec out = arma::join_cols(as<arma::vec>(Opt[0]), as<arma::vec>(Opt[1]));
  
  return out;
}
