#' @title Generate Latent Time Series based on Model
#' @description Create a time series based on a supplied time series model.
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} containing the amount of observations for the time series.
#' @return A \code{lts} object with the structure:
#' \itemize{
#'   \item{"process"} {A \code{vector} that contains names of decomposed and combined processes}
#'   \item{"data"} {A \code{data frame} that contains data for decomposed and combined processes}
#' }
#' @author JJB, Wenchao
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' gen.lts(model)
gen.lts = function(model, N = 1000){
  
  # Do we have a valid model?
  if(!(is(model, "ts.model") || is(model, "gmwm"))){
    stop("model must be created from a ts.model or gmwm object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  if(is(model,"gmwm")){
    model = model$model.hat
  }
  
  # Information Required by GMWM:
  desc = model$desc
  p = length(desc) # p decomposed processes
  obj = model$obj.desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(!model$starting){
    theta = model$theta
    out = .Call('GMWM_gen_lts', PACKAGE = 'GMWM', N, theta, desc, obj)
  }else{
    stop("Need to supply initial values within the ts.model object.")
  }
  
  # (p+1) processes
  comp.desc = c(desc, paste0(desc, collapse = '+'))
  # N by (p+1) data frame
  data = data.frame(matrix(nrow = N, ncol = p+1))
  colnames(data) = comp.desc
  for (i in 1:(p+1)){
    data[[i]] = out[,,i]
  }
  
  out = structure(list(process = comp.desc,
                       data = data), class = 'lts' )
  invisible(out)
}

#plot.lts = function()
#autoplot.lts = function()
#demo.lts = function()