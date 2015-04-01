#' @title Generate Time Series based on Model
#' @description Create a time series based on a supplied time series model.
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param N An \code{interger} containing the amount of observations for the time series.
#' @return A \code{vec} that contains combined time series.
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' What is NOT active:
#' 1. Simulating an ARMA time series
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' gen.ts(model)
gen.ts = function(model, N = 1000){
  
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$obj.desc
  
  obj = model$obj
      
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }

  if(model$adv){
    theta = model$theta
    out = .Call('GMWM_gen_model', PACKAGE = 'GMWM', N, theta, desc, obj)
  }else{
    stop("Need to supply initial values within the ts.model object.")
  }

  out
}