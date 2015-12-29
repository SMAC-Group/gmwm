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

#include "objective_functions.h"

// Uses the theoretical_wv function.
#include "process_to_wv.h"

// Uses the transform / untransform methods
#include "transform_data.h"


// Used Yannick's flattening technique on guessed starting values...
double objFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  // Untransform and find the theoretical wv.
  arma::vec wv_theo = theoretical_wv(untransform_values(theta, desc, objdesc, model_type),
                                     desc, objdesc, tau);
  
  // Yannick's Circle Idea
  arma::vec standardized = 1-wv_theo/wv_empir;
  
  // Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

// Main objective function used by the program
double objFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau){
  
  // Untransform and find the theoretical wv.
  arma::vec wv_theo = theoretical_wv(untransform_values(theta, desc, objdesc, model_type),
                                     desc, objdesc, tau);
  
  // Compute Difference
	arma::vec dif = wv_theo - wv_empir;
	
	// Compute quandratic form... Scalar!!!!!!
	return arma::as_scalar(trans(dif)*omega*dif);
}

//' @title Retrieve GMWM starting value from Yannick's objective function
//' @description Obtains the GMWM starting value given by Yannick's objective function optimization
//' @template tsobj_cpp
//' @param model_type A \code{string} containing the model type. Either 'imu' or 'ssm'
//' @param wv_empir A \code{vec} containing the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales of 2^(1:J), where J is the number of scales created by the decomposition.
//' @return A \code{double} that is the value of the Objective function under Yannick's starting algorithm
//' @keywords internal
// [[Rcpp::export]]
double getObjFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

  return objFunStarting(transformed_theta, desc, objdesc, model_type, wv_empir, tau);
}

//' @title Retrieve GMWM starting value from Yannick's objective function
//' @description Obtains the GMWM starting value given by Yannick's objective function optimization
//' @template tsobj_cpp
//' @param model_type A \code{string} containing the model type. Either 'imu' or 'ssm'
//' @param omega A \code{mat} that is the inverse of the diagonal of the V matrix.
//' @param wv_empir A \code{vec} containing the empirical wavelet variance.
//' @param tau A \code{vec} that contains the scales of 2^(1:J), where J is the number of scales created by the decomposition.
//' @return A \code{double} that is the value of the Objective function under Yannick's starting algorithm
//' @keywords internal
// [[Rcpp::export]]
double getObjFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau){
  
    arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

    return objFun(transformed_theta, desc, objdesc, model_type, omega, wv_empir, tau);
}

/// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                          const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(Rcpp::_["par"] = theta,
                        Rcpp::_["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        Rcpp::_["method"] = "CG",
                        Rcpp::_["desc"] = desc,
                        Rcpp::_["objdesc"] = objdesc,
                        Rcpp::_["model_type"] = model_type,
                        Rcpp::_["wv_empir"] = wv_empir,
                        Rcpp::_["tau"] = tau);
   
   arma::vec out = Rcpp::as<arma::vec>(Opt[0]);
   
   return out;
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, 
                     const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                     const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(Rcpp::_["par"] = theta,
                        Rcpp::_["fn"]  = Rcpp::InternalFunction(&objFun),
                        Rcpp::_["method"] = "CG",
                        Rcpp::_["desc"] = desc,
                        Rcpp::_["objdesc"] = objdesc,
                        Rcpp::_["model_type"] = model_type,
                        Rcpp::_["omega"] = omega,
                        Rcpp::_["wv_empir"] = wv_empir,
                        Rcpp::_["tau"] = tau);
   
   arma::vec out = Rcpp::as<arma::vec>(Opt[0]); 
   
   return out;
}
