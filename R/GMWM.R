#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param data A \code{matrix} or \code{data.frame} object with only column (e.g. \eqn{N \times 1}{ N x 1 })
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU or SSM
#' @param fast A \code{boolean} indicating the type of covariance matrix solver. TRUE means bootstrap and FALSE means use a fast fourier transform (fft)
#' @param augmented A \code{boolean} indicating whether to add additional moments (e.g. mean for drift and variance for all other components).
#' @param p A \code{double} between 0 and 1 that correspondings to the \eqn{\frac{\alpha}{2}}{alpha/2} value for the wavelet confidence intervals.
#' @param robust A \code{boolean} indicating whether to use the robust computation (TRUE) or not (FALSE).
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency.
#' @param B An \code{integer} to sample the space for IMU and SSM models to ensure optimal identitability.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' What is NOT active:
#' 1. Fast Computation of V
#' 2. Augmented GMWM (additional moments)
#' 3. Guessing ARMA models is NOT supported. You must supply specific parameters (see example below)
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' data = gen_ar1(n, phi=.99, sigma2 = 0.01) + gen_wn(n, sigma2=1)
#' 
#' # Models can contain specific parameters e.g.
#' adv.model = gmwm(AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 0.01),
#'                             data)
#' 
#' # Or we can guess the parameters:
#' guided.model = gmwm(AR1() + WN(), data) 
#' 
#' # Want to try different models? 
#' guided.ar1 = gmwm(AR1(), data)
#' 
#' # Faster:
#' guided.ar1.wn.prev = update(guided.ar1, AR1()+WN())
#' 
#' # OR 
#' 
#' # Create new GMWM object. 
#' # Note this is SLOWER since the Covariance Matrix is recalculated.
#' guided.ar1.wn.new = gmwm(AR1()+WN(), data)
#'  
#' # ARMA case
#' set.seed(1336)
#' data = arima.sim(n = 200, 
#'               list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
#'               sd = sqrt(0.1796))
#' guided.arma = gmwm(ARMA(2,2), data, model.type="ssm")
#' adv.arma = gmwm(ARMA(ar=c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), sigma2=0.1796),
#'                 data, model.type="ssm")
gmwm = function(model, data, model.type="imu", fast=TRUE, augmented=FALSE, p = 0.025, robust=FALSE, eff=0.6, B=1000){
  
  # Are we receiving one column of data?
  if( (class(data) == "data.frame" && ncol(data) > 1) || ( class(data) == "matrix" && ncol(data) > 1 ) ){
    stop("The function can only operate on one column of data")
  }
  
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$obj.desc
  
  obj = model$obj
  
  np = model$plength
  
  N = length(data)
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  # Model type issues
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either IMU or SSM!")
  }
  
  # Verify Scales and Parameter Space
  nlevels =  floor(log2(length(data)))
  scales = .Call('GMWM_scales_cpp', PACKAGE = 'GMWM', nlevels)
  
  if(np > length(scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
  }
  
  if(robust){
    np = np+1
    if(np > length(scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }
  }
  
  # Decompose the Series
  modwt.decomp = .Call('GMWM_modwt_cpp', PACKAGE = 'GMWM', data, filter_name = "haar", nlevels, boundary="periodic")
  
  # Obtain the wave variance and chi-square confidence interval
  wv.all = .Call('GMWM_wvar_cpp', PACKAGE = 'GMWM', modwt.decomp, robust, eff, p, "eta3", "haar")
  
  wv.empir = wv.all[,1]
  wv.ci.low = wv.all[,2]
  wv.ci.high = wv.all[,3]
  
  # Create the V Matrix (this is very time consuming...)
  Vout = .Call('GMWM_compute_cov_cpp', PACKAGE = 'GMWM', modwt.decomp, nlevels, "diag", robust, eff)
  
  # Assign the appropriate V (it computes both robust and non-robust)
  if(!robust){
    V = Vout[[1]]
  }else{
    V = Vout[[2]]
    
    ind = 1:np
    
    temp.scales = scales
    temp.V = V
    temp.wv.empir = wv.empir
    
    scales = scales[ind]
    wv.empir = wv.empir[ind]
    V = V[ind,ind]
    
    np = np - 1
  }
  
  # Needed if model contains a drift. 
  expect.diff = .Call('GMWM_mean_diff', PACKAGE = 'GMWM', data)
  
  if(!model$adv){
    theta = .Call('GMWM_guess_initial', PACKAGE = 'GMWM', desc, obj, model.type, np, expect.diff, N, wv.empir, scales, B)
  }else{
    theta = model$theta
  }
  
  out = .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, desc, obj, model.type, V, wv.empir, scales)

  if(robust){
    scales = temp.scales
    V = temp.V
    wv.empir = temp.wv.empir
  }
  
  theo = theoretical_wv(out, desc, obj, scales)
  
  colnames(out) = model$desc
  
    
  out = structure(list(estimate = as.vector(out),
                       wv.empir = wv.empir, 
                       ci_low = wv.ci.low, 
                       ci_high = wv.ci.high, 
                       scales = scales, 
                       theo = theo, 
                       V = V,
                       robust = robust,
                       eff = eff,
                       model.type = model.type,
                       B = B,
                       expect.diff = expect.diff,
                       N = N), class = "gmwm")
  invisible(out)
}


update.gmwm = function(object, model, ...){
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$obj.desc
  
  obj = model$obj
  
  np = model$plength
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }

  wv.empir = object$wv.empir  
  wv.ci.low = object$ci_low
  wv.ci.high = object$ci_high
  V = object$V
  scales = object$scales
  robust = object$robust
  eff = object$eff
  model.type = object$model.type
  B = object$B
  expect.diff = object$expect.diff
  N = object$N
  
  if(np > length(scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
  }
  
  if(robust){
    np = np+1
    if(np > length(scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need one additional scale since robust requires the amount of parameters + 1 to estimate.")
    }
    ind = 1:np
    
    temp.scales = scales
    temp.V = V
    temp.wv.empir = wv.empir
    
    scales = scales[ind]
    wv.empir = wv.empir[ind]
    V = V[ind,ind]
  }
  
  if(!model$adv){
    theta = .Call('GMWM_guess_initial', PACKAGE = 'GMWM', desc, obj, model.type, np, expect.diff, N, wv.empir, scales, B)
  }else{
    theta = model$theta
  }
  
  out = .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, desc, obj, model.type, V, wv.empir, scales)
  
  if(robust){
    scales = temp.scales
    V = temp.V
    wv.empir = temp.wv.empir
  }
  
  theo = theoretical_wv(out, desc, obj, scales)
  
  colnames(out) = model$desc
  
  out = structure(list(estimate = as.vector(out),
                       wv.empir = wv.empir, 
                       ci_low = wv.ci.low, 
                       ci_high = wv.ci.high, 
                       scales = scales, 
                       theo = theo, 
                       V = V,
                       robust = robust,
                       eff = eff,
                       model.type = model.type,
                       B=B), class = "gmwm")
  invisible(out)
}


#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models. Not used if \code{type="arma"}
#' @param wvcov A \code{wvcov} object
#' @param signal A \code{vec} that is the time series being studied
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU, ARMA, or SSM
#' @param B An \code{integer} to sample the space for IMU and SSM models to ensure AR1 identitability.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save.guided.ar1 = guided.gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' save.guided.2ar1 = guided.gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#'  
#' # ARMA case
#' set.seed(1336)
#' x = arima.sim(n = 200, 
#'               list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
#'               sd = sqrt(0.1796))
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save.guided.arma = guided.gmwm(ARMA(2,2),wvcov=out, signal=x, model.type="ssm")
guided.gmwm = function(model, wvcov, signal, model.type="imu", B = 1000){
  
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1, DR, RW, QN, WN, ARMA). ")
  }
  
  desc = model$obj.desc
  
  obj = model$obj
  
  np = model$plength
  
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either IMU or SSM!")
  }
  
  if(length(desc) == 0 ){
    stop("desc must contain a list of components (e.g. AR1, DR, RW, QN, WN).")
  }
  
  if(np > length(wvcov$scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
  }
  
  if(wvcov$robust){
    np = np+1
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }
    ind = 1:np
    
    temp.scales = wvcov$scales
    
    wvcov$scales = wvcov$scales[ind]
    wvcov$wv.empir = wvcov$wv.empir[ind]
    wvcov$V = wvcov$V[ind,ind]
    wvcov$ci_low = wvcov$ci_low[ind]
    wvcov$ci_high = wvcov$ci_high[ind]
  }
  
  
  expect.diff = .Call('GMWM_mean_diff', PACKAGE = 'GMWM', data)
  
  theta = .Call('GMWM_guess_initial', PACKAGE = 'GMWM', desc, obj, model.type, np, expect.diff, length(signal), wv.empir, scales, B)
  
  out = .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, desc, obj, model.type, wvcov$V, wvcov$wv.empir, wvcov$scales)
  
  if(wvcov$robust){
    wvcov$scales = temp.scales
  }
  
  theo = theoretical_wv(out, desc, obj, wvcov$scales)
  
  colnames(out) = model$desc
  
  
  out = structure(list(estimate = as.vector(out),
                       wv.empir = wvcov$wv.empir, 
                       ci_low = wvcov$ci_low, 
                       ci_high = wvcov$ci_high, 
                       scales = wvcov$scales, 
                       theo = theo, 
                       robust = wvcov$robust,
                       eff = wvcov$eff,
                       type = model.type), class = "gmwm")
  invisible(out)
}

#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models. Not used if \code{type="arma"}
#' @param wvcov A \code{wvcov} object
#' @param signal A \code{vec} that is the time series being studied
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU, ARMA, or SSM
#' @param B An \code{integer} to sample the space for IMU and SSM models to ensure AR1 identitability.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save.adv.gmwm.ar1 = adv.gmwm(AR1(phi=.1,sigma2=1), wvcov=out, signal=x, model.type="imu")
#' save.adv.gmwm.2ar1 = adv.gmwm(AR1(phi=.1,sigma2=1)+AR1(phi=0.95,sigma2=.1), wvcov=out, signal=x, model.type="imu")
#'  
#' # ARMA case
#' set.seed(1336)
#' x = arima.sim(n = 200, 
#'               list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
#'               sd = sqrt(0.1796))
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save = adv.gmwm(ARMA(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), sigma2 = 0.1796), wvcov=out, signal=x, model.type="ssm")
adv.gmwm = function(model, wvcov, signal, model.type="imu", B = 1000){
  
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1, DR, RW, QN, WN, ARMA). ")
  }
  
  if(!model$adv){
    stop("The model submitted did NOT contain user-specified values for each term. Please recreate the model with specific coefficients.")
  }
  
  theta = model$theta
  
  desc = model$obj.desc
  
  obj = model$obj
  
  np = model$plength
  
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either IMU or SSM!")
  }
  
  if(length(desc) == 0 ){
    stop("desc must contain a list of components (e.g. AR1, DR, RW, QN, WN).")
  }
  
  if(np > length(wvcov$scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
  }
  
  if(wvcov$robust){
    np = np+1
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }
    ind = 1:np
    
    temp.scales = wvcov$scales
    
    wvcov$scales = wvcov$scales[ind]
    wvcov$wv.empir = wvcov$wv.empir[ind]
    wvcov$V = wvcov$V[ind,ind]
    wvcov$ci_low = wvcov$ci_low[ind]
    wvcov$ci_high = wvcov$ci_high[ind]
  }
  
  out = .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, desc, obj, model.type, wvcov$V, wvcov$wv.empir, wvcov$scales)
  
  if(wvcov$robust){
    wvcov$scales = temp.scales
  }
  
  theo = theoretical_wv(out, desc, obj, wvcov$scales)
  
  colnames(out) = model$desc
  
  
  out = structure(list(estimate = as.vector(out),
                       wv.empir = wvcov$wv.empir, 
                       ci_low = wvcov$ci_low, 
                       ci_high = wvcov$ci_high, 
                       scales = wvcov$scales, 
                       theo = theo, 
                       robust = wvcov$robust,
                       eff = wvcov$eff,
                       type = model.type), class = "gmwm")
  invisible(out)
}

#' @title Wrapper to Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method plot gmwm
#' @param x A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' mod = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' plot(mod)
plot.gmwm = function(x, ...){
  autoplot(x)
}

#' @title Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method autoplot gmwm
#' @param object A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' mod = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' autoplot(mod)
autoplot.gmwm = function(object, ...){
  
  low=high=emp=theo=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci_low,
                  high = object$ci_high,
                  scale = object$scales,
                  theo = object$theo)
  CI = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3")) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Haar Wavelet Variance Representation") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta),")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0)) 
  
  CI
}

#' @title Compare GMWM Model Fits on Same Graph
#' @description Creates a single graph that contains two GMWM models plotted against each other.
#' @method autoplot comp
#' @param object A \code{data.frame} containing both sets of GMWM object data.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing one graphs with two GMWM models plotted against each other.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(1335)
#' n = 200
#' x = gen_ar1(n, phi = .1, sigma2 = 1) + gen_ar1(n, phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
autoplot.comp = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#EF8A1C","LINE4"="#6E8BBF")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci_low,
                  high = object$ci_high,
                  scale = object$scales,
                  theo = object$theo1,
                  model2 = object$theo2)
  CI = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3"), size = 1) + 
    geom_line(aes(y = model2, colour = "LINE4"), size = 1) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    geom_point(aes(y = model2, colour = "LINE4"), size = 4, shape = 0) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Haar Wavelet Variance Representation") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Model 1: ", nu,"(",hat(theta)[1],")")),expression(paste("Model 2: ", nu,"(",hat(theta)[2],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1,1), colour = c("#000000","#999999","#EF8A1C","#6E8BBF")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0)) 
  
  CI
}


#' @title Compare GMWM Model Fit on Split Graphs
#' @description Creates GMWM model fits of two models on split graphs within the same panel.
#' @method autoplot compSplit
#' @param object A \code{data.frame} containing both sets of GMWM object data.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing two graphs of the wavelet variance.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(1355)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = TRUE)}
autoplot.compSplit = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci_low,
                  high = object$ci_high,
                  scale = object$scales,
                  theo = object$theo1,
                  model2 = object$theo2)
  
  CI1 = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3")) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Model 1") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta)[1],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0))
  
  CI2 = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = model2, colour = "LINE3")) + 
    geom_point(aes(y = model2, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Model 2") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta)[2],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0))
  
  multiplot(CI1, CI2, cols=2)	
}

#' @title Graphically Compare GMWM Model Fit Between Two Models
#' @description Creates GMWM model fits of two models on split graphs within the same panel.
#' @usage compare.models(GMWM1, GMWM2, split = FALSE)
#' @param GMWM1 A \code{gmwm} object
#' @param GMWM2 A \code{gmwm} object
#' @param split A \code{boolean} indicating true or false to place model fits on different or the same graphs.
#' @return A ggplot2 panel containing two graphs of the wavelet variance.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(8836)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
compare.models = function(GMWM1, GMWM2, split = FALSE){
  x = data.frame(wv.empir = GMWM1$wv.empir, ci_low = GMWM1$ci_low, 
                 ci_high = GMWM1$ci_high, scales = GMWM1$scales, theo1 = GMWM1$theo, theo2 = GMWM2$theo) 
  if (split == TRUE){
    class(x) = "compSplit"
  }else{
    class(x) = "comp"
  }
  autoplot(x)
}
