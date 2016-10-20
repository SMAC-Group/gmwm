/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
 
#include "ts_model_cpp.h"

// Goal of file is to be able to recreate a string of values given.
// Not supported: ARMA(P,Q) model
// Support does exist for the ARMA(1,1) model


//' @title Transform AR1 to GM
//' @description 
//' Takes AR1 values and transforms them to GM
//' @param theta A \code{vec} that contains AR1 values.
//' @param freq  A \code{double} indicating the frequency of the data.
//' @return A \code{vec} containing GM values.
//' @details
//' The function takes a vector of AR1 values \eqn{\phi}{phi} and \eqn{\sigma ^2}{sigma ^2}
//' and transforms them to GM values \eqn{\beta}{beta} and \eqn{\sigma ^2_{gm}}{sigma ^2[gm]}
//' using the formulas:
//' \eqn{\beta  =  - \frac{{\ln \left( \phi  \right)}}{{\Delta t}}}{beta = -ln(phi)/delta_t}
//' \eqn{\sigma _{gm}^2 = \frac{{{\sigma ^2}}}{{1 - {\phi ^2}}} }{sigma^2[gm] = sigma^2/(1-phi^2)}
//' @keywords internal
//' @author JJB
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
//' @examples
//' ar1_to_gm(c(0.3,1,0.6,.3), 2)
// [[Rcpp::export]]
arma::vec ar1_to_gm(arma::vec theta, double freq){
  unsigned int n = theta.n_elem;
  if(theta.n_elem %2 != 0){Rcpp::stop("Bad Theta Vector");}
  
  // Convert to delta.t
  freq = 1.0/freq;
  
  for(unsigned int i = 0; i < int(double(n)/2.0); i++){
    double phi = theta(2*i);
    double sigma2 = theta(2*i+1);
    theta(2*i) = -1*log(phi)/freq;         // beta 
    theta(2*i+1) = sigma2 / (1.0-phi*phi); // sigma2_gm
  }
  
  return theta;
}


//' @title Transform GM to AR1
//' @description Takes GM values and transforms them to AR1
//' @param theta A \code{vec} that contains AR1 values.
//' @param freq A \code{double} indicating the frequency of the data.
//' @return A \code{vec} containing GM values.
//' @keywords internal
//' @author JJB
//' The function takes a vector of GM values \eqn{\beta}{beta} and \eqn{\sigma ^2_{gm}}{sigma ^2[gm]}
//' and transforms them to AR1 values \eqn{\phi}{phi} and \eqn{\sigma ^2}{sigma ^2}
//' using the formulas:
//' \eqn{\phi  = \exp \left( { - \beta \Delta t} \right)}{phi = exp(-beta * delta[t])}
//' \eqn{{\sigma ^2} = \sigma _{gm}^2\left( {1 - \exp \left( { - 2\beta \Delta t} \right)} \right)}{sigma^2 = sigma^2[gm]*(1-exp(-2*beta*delta[t]))}
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
//' @examples
//' gm_to_ar1(c(0.3,1,0.6,.3), 2)
// [[Rcpp::export]]
arma::vec gm_to_ar1(arma::vec theta, double freq){
  unsigned int n = theta.n_elem;
  if(theta.n_elem %2 != 0){Rcpp::stop("Bad Theta Vector");}
  
  // Convert to delta.t
  freq = 1.0/freq;
  
  for(unsigned int i = 0; i < int(double(n)/2.0); i++){
    double beta = theta(2*i);
    double sigma2_gm = theta(2*i+1);
    theta(2*i) = exp(-1*beta*freq);           // phi 
    theta(2*i+1) = sigma2_gm*(1-exp(-2*beta*freq)); // sigma2
  }
  
  return theta;
}


 
//' @title Generate the ts model object description
//' @description Creates the ts.model's obj.desc value
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{field<vec>} that contains the object description of each process.
//' @details
//' This function currently does NOT support ARMA(P,Q) models. 
//' That is, there is no support for ARMA(P,Q), AR(P), or MA(Q).
//' There is support for ARMA11, AR1, MA1, GM, WN, DR, QN, and RW.
//' @keywords internal
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
// [[Rcpp::export]]
arma::field<arma::vec> model_objdesc(std::vector<std::string> desc){
  unsigned int n = desc.size();
  arma::field<arma::vec> objdesc(n);
  
  arma::vec ar1 = arma::ones<arma::vec>(2);
  arma::vec arma11 = arma::ones<arma::vec>(3);
  arma::vec others = arma::ones<arma::vec>(1);
  
  for(unsigned int i = 0; i < n; i++){
    std::string element_type = desc[i];
    if(element_type == "AR1" || element_type == "GM" || element_type == "MA1"){
      objdesc(i) = ar1;
    }else if(element_type != "ARMA11"){
      objdesc(i) = others;
    }else{
      objdesc(i) = arma11;
    }
  }
  
  return objdesc;
}


//' @title Generate the ts model object's theta vector
//' @description Creates the ts.model's theta vector
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{vec} with values initialized at 0 that span the space of parameters to be estimated.
//' @details
//' This function currently does NOT support ARMA(P,Q) models. 
//' That is, there is no support for ARMA(P,Q), AR(P), or MA(Q).
//' There is support for ARMA11, AR1, MA1, GM, WN, DR, QN, and RW.
//' @keywords internal
//' @backref src/ts_model_cpp.cpp
//' @backref src/ts_model_cpp.h
// [[Rcpp::export]]
arma::vec model_theta(std::vector<std::string> desc){
  unsigned int n = desc.size();
  
  unsigned int m = 0;
  for(unsigned int i = 0; i < n; i++){
    std::string element_type = desc[i];
    if(element_type == "AR1" || element_type == "GM" || element_type == "MA1"){
      m += 2;
    }else if(element_type != "ARMA11"){
      m++;
    }else{
      m += 3;
    }
  }
  
  return arma::zeros<arma::vec>(m);
}

//' @title Generate the ts model object's process desc
//' @description Creates the ts.model's process desc
//' @param desc A \code{vector<string>} that contains a list of the strings of each process.
//' @return A \code{vector<string>} with a list of descriptive values to label the estimate matrix with
//' @details
//' This function currently does NOT support ARMA(P,Q) models. 
//' That is, there is no support for ARMA(P,Q), AR(P), or MA(Q).
//' There is support for ARMA11, AR1, MA1, GM, WN, DR, QN, and RW.
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
    }else if(element_type == "GM"){
      proc_desc.push_back("BETA");
      proc_desc.push_back("SIGMA2_GM");
    }else if(element_type == "MA1"){
      proc_desc.push_back("MA1");
      proc_desc.push_back("SIGMA2");
    }else if(element_type != "ARMA11"){
      proc_desc.push_back(element_type);
    }else{
      proc_desc.push_back("AR1");
      proc_desc.push_back("MA1");
      proc_desc.push_back("SIGMA2");
    }
  }
  
  return proc_desc;
}
