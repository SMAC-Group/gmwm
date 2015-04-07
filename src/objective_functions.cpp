#include <RcppArmadillo.h>

#include "objective_functions.h"
#include "process_to_wv.h"
#include "transform_data.h"

using namespace Rcpp;


// Used Yannick's flattening technique on guessed starting values...
double objFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec untransformed_theta = untransform_values(theta, desc, objdesc, model_type);
  arma::vec wv_theo = theoretical_wv(untransformed_theta, desc, objdesc, tau);
  arma::vec standardized = 1-wv_theo/wv_empir;
  // Compute quandratic form
	return arma::as_scalar(trans(standardized)*(standardized));
}

// Main objective function used by the program
double objFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega,const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec untransformed_theta = untransform_values(theta, desc, objdesc, model_type);

  arma::vec wv_theo = theoretical_wv(untransformed_theta, desc, objdesc, tau);

  // Compute quandratic form
	arma::vec dif = wv_theo - wv_empir;
	return arma::as_scalar(trans(dif)*omega*dif);
}

// [[Rcpp::export]]
double getObjFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau){
  
  arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

  return objFunStarting(transformed_theta, desc, objdesc, model_type, wv_empir, tau);
}

// [[Rcpp::export]]
double getObjFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& V,const arma::vec& wv_empir, const arma::vec& tau){
  
    arma::vec transformed_theta = transform_values(theta, desc, objdesc, model_type);

    return objFun(transformed_theta, desc, objdesc, model_type, arma::inv(V), wv_empir, tau);
}

/// [[Rcpp::export]]
arma::vec Rcpp_OptimStart(const arma::vec&  theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                          const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFunStarting),
                        _["method"] = "CG",
                        _["desc"] = desc,
                        _["objdesc"] = objdesc,
                        _["model_type"] = model_type,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}

/// [[Rcpp::export]]
arma::vec Rcpp_Optim(const arma::vec&  theta, 
                     const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                     const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau){
   Rcpp::Environment stats("package:stats"); 
   Rcpp::Function optim = stats["optim"];    

   Rcpp::List Opt=optim(_["par"] = theta,
                        _["fn"]  = Rcpp::InternalFunction(&objFun),
                        _["method"] = "CG",
                        _["desc"] = desc,
                        _["objdesc"] = objdesc,
                        _["model_type"] = model_type,
                        _["omega"] = omega,
                        _["wv_empir"] = wv_empir,
                        _["tau"] = tau);
   
   arma::vec out = as<arma::vec>(Opt[0]);
   
   return out;
}
