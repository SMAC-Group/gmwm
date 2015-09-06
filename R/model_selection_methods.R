#' @title Create Combination Matrix
#' @description Generates a matrix containing all the combination listings
#' @param n The number of variables.
#' @details Port expand.grid to C++ at a later time...
#' @return Returns a binary matrix (e.g. 1 or 0) entries
comb.mat = function(n){
  c = rep(list(1:0), n)
  expand.grid(c)
}

#' @title Automatically select appropriate model for IMU
#' @description Runs through a model selection algorithm to determine the best model
#' @param model A \code{ts.model} object that is the largest model to be tested.
#' @param data A \code{vector}, \code{matrix}, \code{data.frame}, or \code{imu} object with either 1, 3, or 6 columns. 
#' @param bootstrap A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param alpha A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust A \code{boolean} that indicates whether to use robust estimation.
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B A \code{integer} that contains the amount of bootstrap replications
#' @param G A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @return A \code{auto.imu} object.
auto.imu = function(data, model = 3*AR1()+WN()+RW()+QN()+DR(), bootstrap = F, alpha = 0.05, robust = F, eff = 0.6, B = 20, G = 1000){
  
  data = as.matrix(data)
  full.str = model$desc
  m = as.matrix(comb.mat(length(full.str)))
  m = m[-nrow(m),]
  
  out = .Call('GMWM_auto_imu', PACKAGE = 'GMWM', data, combs=m, full_model=full.str, alpha, compute_v = "fast", model_type = "sensor", K=1, H=B, G, robust, eff, bootstrap)
  
}