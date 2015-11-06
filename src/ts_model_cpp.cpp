/* Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify it
 * under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
 * as the LICENSE file.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
 * (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.
 * 
 */

#include <RcppArmadillo.h>
using namespace Rcpp;
 
#include "ts_model_cpp.h"

// Goal of file is to be able to recreate a string of values given.
// Not supported: ARMA
 
//' @title Generate the ts model object description
//' @description Creates the ts.model's obj.desc value
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{field<vec>} that contains the object description of each process.
//' @details
//' This function currently does NOT support ARMA models. 
//' That is, there is no support for ARMA, AR, or MA.
//' There is support for AR1, WN, DR, QN, and RW.
//' @keywords internal
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
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


//' @title Generate the ts model object's theta vector
//' @description Creates the ts.model's theta vector
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{vec} with values initialized at 0 that span the space of parameters to be estimated.
//' @details
//' This function currently does NOT support ARMA models. 
//' That is, there is no support for ARMA, AR, or MA.
//' There is support for AR1, WN, DR, QN, and RW.
//' @keywords internal
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
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

//' @title Generate the ts model object's process desc
//' @description Creates the ts.model's process desc
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{vector<string>} with a list of descriptive values to label the estimate matrix with
//' @details
//' This function currently does NOT support ARMA models. 
//' That is, there is no support for ARMA, AR, or MA.
//' There is support for AR1, WN, DR, QN, and RW.
//' @keywords internal
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
// [[Rcpp::export]]
std::vector<std::string> model_process_desc(std::vector<std::string> desc){
  unsigned int n = desc.size();

  std::vector<std::string> proc_desc;
  
  for(unsigned int i = 0; i < n; i++){
    std::string element_type = desc[i];
    if(element_type == "AR1"){
      proc_desc.push_back("AR1");
      proc_desc.push_back("SIGMA2");
    }else{
      proc_desc.push_back(element_type);
    }
  }
  
  return proc_desc;
}
