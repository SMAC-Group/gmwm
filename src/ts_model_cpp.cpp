#include <RcppArmadillo.h>
using namespace Rcpp;
 
#include "ts_model_cpp.h"

// Goal of file is to be able to recreate a string of values given.
// Not supported: ARMA
 
 
 // [[Rcpp::export]]
arma::field<arma::vec> model_objdesc(std::vector<std::string> desc){
  unsigned int n = desc.size();
  arma::field<arma::vec> objdesc(n);
  
  arma::vec ar1 = arma::ones<arma::vec>(2);
  arma::vec others = arma::ones<arma::vec>(1);
  
  for(unsigned int i = 0; i < n; i++){
    std::string element_type = desc[i];
    if(element_type == "AR1"){
      objdesc(i) = ar1;
    }else{
      objdesc(i) = others;
    }
  }
  
  return objdesc;
}
 
 
// [[Rcpp::export]]
arma::vec model_theta(std::vector<std::string> desc){
  unsigned int n = desc.size();
  
  unsigned int m = 0;
  for(unsigned int i = 0; i < n; i++){
    std::string element_type = desc[i];
    if(element_type == "AR1"){
      m += 2;
    }else{
      m++;
    }
  }
  
  return arma::zeros<arma::vec>(m);
}
