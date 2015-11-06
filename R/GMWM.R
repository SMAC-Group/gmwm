# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.


#' @title GMWM for Sensors, ARMA, SSM, and Robust
#' @description GMM object
#' @export
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param data A \code{matrix} or \code{data.frame} object with only column (e.g. \eqn{N \times 1}{ N x 1 }), or a \code{lts} object, or a \code{gts} object. 
#' @param model.type A \code{string} containing the type of GMWM needed e.g. sensor or SSM
#' @param compute.v A \code{string} indicating the type of covariance matrix solver. "fast", "bootstrap", "asymp.diag", "asymp.comp", "fft"
#' @param augmented A \code{boolean} indicating whether to add additional moments (e.g. mean for drift and variance for all other components).
#' @param alpha A \code{double} between 0 and 1 that correspondings to the \eqn{\frac{\alpha}{2}}{alpha/2} value for the wavelet confidence intervals.
#' @param robust A \code{boolean} indicating whether to use the robust computation (TRUE) or not (FALSE).
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency.
#' @param G An \code{integer} to sample the space for sensor and SSM models to ensure optimal identitability.
#' @param K An \code{integer} that controls how many times the bootstrapping procedure will be initiated.
#' @param H An \code{integer} that indicates how many different samples the bootstrap will be collect.
#' @param seed A \code{integer} that controls the reproducibility of the auto model selection phase.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' What is NOT active:
#' 1. Augmented GMWM (additional moments)
#' 
#' The V matrix is calculated by:
#' \eqn{diag\left[ {{{\left( {Hi - Lo} \right)}^2}} \right]}{diag[(Hi-Lo)^2]}.
#' 
#' The function is implemented in the following manner:
#' 1. Calculate MODWT of data with levels = floor(log2(data))
#' 2. Apply the brick.wall of the MODWT (e.g. remove boundary values)
#' 3. Compute the empirical wavelet variance (WV Empirical).
#' 4. Obtain the V matrix by squaring the difference of the WV Empirical's Chi-squared confidence interval (hi - lo)^2
#' 5. Optimize the values to obtain \eqn{\hat{\theta}}{theta^hat}
#' 6. If FAST = TRUE, return these results. Else, continue.
#' 
#'Loop  k = 1 to K
#' Loop h = 1 to H
#' 7. Simulate xt under \eqn{F_{\hat{\theta}}}{F_theta^hat}
#' 8. Compute WV Empirical
#' END
#' 9. Calculate the covariance matrix
#' 10. Optimize the values to obtain \eqn{\hat{\theta}}{theta^hat}
#'END
#' 11. Return optimized values.
#' 
#' 
#' The function estimates a variety of time series models. If type = "imu" or "SSM", then
#' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
#' options are:
#' \itemize{
#'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
#'   \item{"ARMA"}{an autoregressive moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
#'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
#'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
#'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
#'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
#' }
#' If type = "ARMA", the function takes condition least squares as starting values; if type = "imu" or type = "SSM" then
#' starting values pass through an initial bootstrap and pseudo-optimization before being passed to the GMWM optimization.
#' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
#' 
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' data = gen.gts(AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 1), n)
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
#' #guided.arma = gmwm(ARMA(2,2), data, model.type="ssm")
#' adv.arma = gmwm(ARMA(ar=c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), sigma2=0.1796),
#'                 data, model.type="ssm")
gmwm = function(model, data, model.type="ssm", compute.v="auto", augmented=FALSE, robust=FALSE, eff=0.6, alpha = 0.05, seed = 1337, G = NULL, K = 1, H = 100){
  
  #check lts
  if(is(data,'lts')){
    warning('lts object is detected. This function can only operate on the combined process.')
    data = data$data[ ,ncol(data$data)]
  }
  
  #check gts
  if(is(data,'gts')){
    data = data$data[,1]
  }
  
  # Are we receiving one column of data?
  if( (class(data) == "data.frame" && ncol(data) > 1) || ( class(data) == "matrix" && ncol(data) > 1 ) ){
    stop("The function can only operate on one column of data")
  }
  
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), AR(p), MA(q), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  np = model$plength
  
  N = length(data)
  
  starting = model$starting
  
  # Input guessing
  if((is.null(G) & starting) || !is.whole(G)){
    if(N > 10000){
      G = 1000
    }else{
      G = 20000
    }
  }else if(!starting){
    G = 0
  }
  
  # For reproducibility
  set.seed(seed)
  
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  # Model type issues
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either sensor or imu!")
  }
  
  # Verify Scales and Parameter Space
  nlevels =  floor(log2(length(data)))
  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)
  
  if(np > length(scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need at least the same number of scales as parameters to estimate.")
  }
  
  if(robust){
    np = np+1
    if(np > length(scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we  at least need the same number of scales as parameters to estimate.")
    }
  }

  
  # Auto setting

  # Compute fast covariance if large sample, otherwise, bootstrap.
  if(compute.v == "auto" || ( compute.v != "fast" && compute.v != "diag")){
    compute.v = "fast"
  }
  

  theta = model$theta

  out = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj, model.type, starting = model$starting,
                                                         p = alpha, compute_v = compute.v, K = K, H = H, G = G,
                                                         robust=robust, eff = eff)
  #colnames(out) = model$desc
  
  estimate = out[[1]]
  rownames(estimate) = model$process.desc
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  
  # Create a new model object.
  model.hat = model
  
  model.hat$starting = F  
  
  model.hat$theta = as.numeric(estimate)
  
  # Release model
  out = structure(list(estimate = estimate,
                       init.guess = init.guess,
                       wv.empir = out[[3]], 
                       ci.low = out[[4]], 
                       ci.high = out[[5]],
                       orgV = out[[7]],
                       V = out[[6]],
                       omega = out[[12]],
                       obj.fun = out[[11]],
                       theo = out[[9]],
                       decomp.theo = out[[10]],
                       scales = scales, 
                       robust = robust,
                       eff = eff,
                       model.type = model.type,
                       compute.v = compute.v,
                       augmented = augmented,
                       alpha = alpha,
                       expect.diff = out[[8]],
                       N = N,
                       G = G,
                       H = H,
                       K = K,
                       model = model,
                       model.hat = model.hat,
                       starting = model$starting,
                       seed = seed), class = "gmwm")
  invisible(out)
}

#' @title Update GMWM object for sensor, ARMA, SSM, and Robust
#' @description GMM object
#' @export
#' @param object A \code{gmwm} object.
#' @param model A \code{ts.model} object containing one of the allowed models
#' @param ... Additional parameters (not used)
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' exact.model = AR1(phi=.99, sigma2 = 0.01) + WN(sigma2=1)
#' data = gen.gts(exact.model)
#' 
#' # Create an initial model that is not accurate
#' bad.model = gmwm(AR1(), data = data)
#' 
#' # Models can contain specific parameters e.g.
#' updated.model = update(bad.model, exact.model)
#' 
#' # Or...
#' updated.model.guided = update(bad.model, AR1()+AR1())
update.gmwm = function(object, model, ...){
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  np = model$plength
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  models.active = count_models(desc)
  
  # Set seed for reproducibility
  
  set.seed(object$seed)
  
  # Identifiability issues
  if(any( models.active[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  # ID error:
  if( sum(models.active) == 1 & models.active["ARMA"] == 1 & model$starting){
    warning("ARMA starting guesses using update.gmwm are NOT based on CSS but an alternative algorithm.")
  }
  
  if(np > length(object$scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need  at least  the same number of scales as parameters to estimate.")
  }
  
  if(object$robust){
    np = np+1
    if(np > length(object$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need one additional scale since robust requires the amount of parameters + 1 to estimate.")
    }
  }
  
  out = .Call('gmwm_gmwm_update_cpp', PACKAGE = 'gmwm',
                  model$theta,
                  desc, obj, 
                  object$model.type, object$N, object$expect.diff, 
                  object$orgV, object$scales, object$wv.empir,
                  model$starting, 
                  object$compute.v, object$K, object$H,
                  object$G, 
                  object$robust, object$eff)

  print(out)
  estimate = out[[1]]
  
  model.hat = model
  
  model.hat$starting = F  
  
  model.hat$theta = as.numeric(estimate)
  
  object$model.hat = model.hat
  
  rownames(estimate) = model$process.desc
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  
  object$estimate = estimate
  object$init.guess = init.guess
  
  object$V = out[[3]]
  object$theo = out[[4]]
  object$decomp.theo = out[[5]]
  
  object$starting = model$starting

  invisible(object)
}


#' @title GMWM for (Robust) sensor
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param data A \code{matrix} or \code{data.frame} object with only column (e.g. \eqn{N \times 1}{ N x 1 }), or a \code{lts} object, or a \code{gts} object. 
#' @param compute.v A \code{string} indicating the type of covariance matrix solver. "fast", "bootstrap", "asymp.diag", "asymp.comp", "fft"
#' @param robust A \code{boolean} indicating whether to use the robust computation (TRUE) or not (FALSE).
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency.
#' @param ... other arguments passed to the main gmwm function
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
gmwm.imu = function(model, data, compute.v = "fast", robust = F, eff = 0.6, ...){
  
  x = gmwm(model = model, 
       data = data, 
       compute.v = compute.v,
       model.type = "imu",
       robust = robust, 
       eff = eff,
       ...
       )
  class(x) = c("gmwm_imu","gmwm")
  
  x
}

#' @title Summary of GMWM object
#' @description Displays summary information about GMWM object
#' @method summary gmwm
#' @export
#' @param object A \code{GMWM} object
#' @param inference A value containing either: NULL (auto), TRUE, or FALSE
#' @param bs.gof A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.gof.p.ci A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.theta.est A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.ci A value containing either: NULL (auto), TRUE, FALSE
#' @param B An \code{int} that indicates how many bootstraps should be performed.
#' @param ... other arguments passed to specific methods
#' @return A \code{summary.gmwm} object with:
#' \itemize{
#'  \item{estimates}{Estimated Theta Values}
#'  \item{testinfo}{Goodness of Fit Information}
#'  \item{inference}{Inference performed? T/F}
#'  \item{bs.gof}{Bootstrap GOF? T/F}
#'  \item{bs.gof.p.ci}{Bootstrap GOF P-Value CI? T/F}
#'  \item{bs.theta.est}{Bootstrap Theta Estimates? T/F}
#'  \item{bs.ci}{Bootstrap CI? T/F}
#'  \item{starting}{Indicates if program supplied initial starting values}
#'  \item{seed}{Seed used during guessing / bootstrapping}
#'  \item{obj.fun}{Value of obj.fun at minimized theta}
#'  \item{N}{Length of Time Series}
#' }
#' @author JJB
#' @examples
#' \dontrun{
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen.gts(AR1(phi=.1, sigma2 = 1) + AR1(phi=0.95, sigma2 = .1),n)
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' summary(mod)
#' }
summary.gmwm = function(object, inference = NULL,  
                        bs.gof = NULL,  bs.gof.p.ci = NULL, 
                        bs.theta.est = NULL, bs.ci = NULL,
                        B = 100, ...){
  
  # Set a different seed to avoid dependency.
  set.seed(object$seed+5)
  
  out = object$estimate
  colnames(out) = c("Estimates") 
  
  N = object$N
  
  # Enable values if small time series.
  auto = if(N > 10000) FALSE else TRUE
  
  # Auto set values
  if(is.null(inference)){
    inference = auto
  }

  if(is.null(bs.gof)){
    bs.gof= if(inference) auto else F
  }
  
  if(is.null(bs.gof.p.ci)){
    bs.gof.p.ci = if(inference) auto else F
  }
  
  if(is.null(bs.theta.est)){
    bs.theta.est = if(inference) auto else F
  }
  
  if(is.null(bs.ci)){
    bs.ci = if(inference) auto else F
  }
  
  if("ARMA" %in% object$model$desc){
    if(bs.ci == FALSE){
      warning(paste0("The numerical derivative of ARMA(p,q), where p > 1 and q > 1, may be inaccurate leading to inappropriate CIs.\n",
              "Consider using the bs.ci = T option on the summary function."))
    }
  }
  
  if(inference){
    mm = .Call('gmwm_get_summary', PACKAGE = 'gmwm',object$estimate,
                                                    object$model$desc, object$model$obj.desc,
                                                    object$model.type, 
                                                    object$wv.empir, object$theo,object$scales,
                                                    object$V, solve(object$orgV), object$obj.fun,
                                                    N, object$alpha,
                                                    object$robust, object$eff,
                                                    inference, F, # fullV is always false. Need same logic updates.
                                                    bs.gof, bs.gof.p.ci, bs.theta.est, bs.ci,
                                                    B)
  }else{
    mm = vector('list',3)
    mm[1:3] = NA
  }
  
  if(inference){
    out.coln = colnames(out)
    out = cbind(out, mm[[1]])
    colnames(out) = c(out.coln, "CI Low", "CI High", "SE")
  }
  
  x = structure(list(estimates=out, 
                     testinfo=mm[[2]],
                     inference = inference, 
                     bs.gof = bs.gof,
                     bs.gof.p.ci = bs.gof.p.ci,
                     bs.theta.est = bs.theta.est, 
                     bs.ci = bs.ci,
                     starting = object$starting,
                     seed = object$seed,
                     obj.fun = object$obj.fun,
                     N = N), class = "summary.gmwm")
    
  x
}

#' @title Print summary.gmwm object
#' @description Displays summary information about GMWM object
#' @method print summary.gmwm
#' @export
#' @keywords internal
#' @param x A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' \dontrun{
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen.gts(AR1(phi=.1, sigma2 = 1) + AR1(phi=0.95, sigma2 = .1),n)
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' summary(mod)
#' }
print.summary.gmwm = function(x, ...){
  
  cat("Model Information: \n")
  print(x$estimates)
  if(x$bs.theta.est){
    cat("\n> The parameter estimates shown are bootstrapped! To use these results, please save the summary object.")
  }
  
  cat("\n* The initial values of the parameters used in the minimization of the GMWM objective function were ", 
      {if(x$starting) paste0("generated by the program underneath seed: ",x$seed,".") else "given by YOU!"},"\n\n")

  cat(paste0("Objective Function: ", round(x$obj.fun,4),"\n\n"))
  
    
  if(x$inference){
    cat(paste0({if(x$bs.gof) "Bootstrapped" else "Asymptotic"}," Goodness of Fit: \n"))
    if(x$bs.gof){
      cat(paste0("Test Statistic: ", round(x$obj.fun,2),"\n",
                 "P-Value: ", round(x$testinfo[1],4)), 
                {if(x$bs.gof.p.ci) paste0(" CI: (", round(x$testinfo[2],4),", ", round(x$testinfo[3],4), ")")})

    }else{
      cat(paste0("Test Statistic: ", round(x$testinfo[1],2),
          " on ",x$testinfo[3]," degrees of freedom\n",
          "The resulting p-value is: ", round(x$testinfo[2],4)))
    }
    cat("\n\n")
  }
  
  if(x$bs.gof || x$bs.theta.est)
   cat(paste0("\nTo replicate the results, use seed: ",x$seed, "\n"))
}


#' @title Predict future points in the time series using the solution of the Generalized Method of Wavelet Moments
#' @description Creates a prediction using the estimated values of GMWM through the ARIMA function within R.
#' @method predict gmwm
#' @export
#' @param object A \code{gmwm} object 
#' @param data.in.gmwm The data SAME EXACT DATA used in the GMWM estimation
#' @param n.ahead Number of observations to guess.
#' @param ... Additional parameters
#' @return A \code{predict.gmwm} object with:
#' \itemize{
#' \item{pred}{Predictions}
#' \item{se}{SE}
#' \item{resid}{Residuals}
#' }
predict.gmwm = function(object, data.in.gmwm, n.ahead = 1, ...){
  
  ts.mod = object$model
  
  if(length(ts.mod$desc) > 1 || ts.mod$desc != "ARMA")
    stop("The predict function only works with stand-alone ARMA models.")
  
  objdesc = ts.mod$obj.desc[[1]]
  
  p = objdesc[1]
  q = objdesc[2]
  
  mod = arima(data.in.gmwm, order = c(p, 0, q),
              method="ML",
              fixed = object$estimate[1:(p+q)],
              transform.pars = F,
              include.mean = F)
  
  pred = predict(mod, n.ahead = n.ahead, newxreg = NULL,
                 se.fit = TRUE, ...)
  
  
  out = structure(list(pred = pred$pred,
                       se = pred$se,
                       resid = mod$residuals)
                  , class = "predict.gmwm")
                          
}

#' @title Wrapper to Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method plot gmwm
#' @export
#' @param x A \code{GMWM} object
#' @param process.decomp A \code{boolean} that indicates whether the decomposed processes should be plotted or not
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @template CommonParams
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB, Wenchao
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' mod = gmwm(AR1(), data=x, model.type="imu")
#' plot(mod)
plot.gmwm = function(x, process.decomp = FALSE, background = 'white', CI = T, transparence = 0.1, bw = F, 
                         CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                         point.size = NULL,point.shape = NULL,
                         title = NULL, title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)),
                         legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ... ){
  
  autoplot.gmwm(x, process.decomp = process.decomp, background = background, CI = CI, transparence = transparence, bw = bw, 
                 CI.color = CI.color, line.type = line.type, line.color = line.color,
                 point.size = point.size, point.shape = point.shape,
                 title = title, title.size= title.size, 
                 axis.label.size = axis.label.size, axis.tick.size = axis.tick.size , 
                 axis.x.label = axis.x.label,
                 axis.y.label = axis.y.label,
                 legend.title = legend.title,  legend.label = legend.label, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                 legend.text.size = legend.text.size)
  
}


#' @title Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method autoplot gmwm
#' @export
#' @keywords internal
#' @param object A \code{GMWM} object
#' @param process.decomp A \code{boolean} that indicates whether the decomposed processes should be plotted or not
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @template CommonParams
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB, Wenchao
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' mod = gmwm(AR1(), data=x, model.type="imu")
#' autoplot(mod)
#' 
#' y = gen.gts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' mod = gmwm(2*AR1(), data = y)
#' autoplot(mod)
autoplot.gmwm = function(object, process.decomp = FALSE, background = 'white', CI = T, transparence = 0.1, bw = F, 
                     CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                     point.size = NULL,point.shape = NULL,
                     title = NULL, title.size= 15, 
                     axis.label.size = 13, axis.tick.size = 11, 
                     axis.x.label = expression(paste("Scale ", tau)),
                     axis.y.label = expression(paste("Wavelet Variance ", nu)),
                     legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                     legend.text.size = 13, ... ){
  
  ## check parameters
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  if(!process.decomp){
    if(CI == T) {numLabel = 3}else {numLabel = 2}
  }else{
    L = length(object$model$desc) + 1 # Find number of latent processes
    if(!bw && (L-1)> 9){warning('Object has more than 9 latent processes, but the palette has only 9 colors')}
    if(CI == T) {numLabel = 2+L}else{numLabel = 1+L}
  }
  
  params = c('line.type', 'line.color', 'point.size','point.shape','legend.label')
  for(i in 1:length(params)){
    one_param = params[i]
    if( !is.null(get(one_param)) && length(get(one_param))!=numLabel){
      warning(paste('Parameter', one_param, 'requires',numLabel,'elements,','but', length(get(one_param)),
                    'is supplied.','Default setting is used.'))
      assign(one_param,NULL)
    }
  }
  
  if(!is.null(legend.label) && anyDuplicated(legend.label) >0){
    warning('Parameter legend.label contains duplicate elements. Add white spaces to each element to avoid this.')
    for(i in 1:length(legend.label) ){
      legend.label[i] = paste0(legend.label[i], paste0(rep(' ',times = i), collapse = ''))
    }
  }
  
  ## call
  if(process.decomp){
    # plot individually
    autoplot.gmwm2(object, background = background, CI = CI, transparence = transparence, bw = bw, 
                   CI.color = CI.color, line.type = line.type, line.color = line.color,
                   point.size = point.size, point.shape = point.shape,
                   title = title, title.size= title.size, 
                   axis.label.size = axis.label.size, axis.tick.size = axis.tick.size , 
                   axis.x.label = axis.x.label,
                   axis.y.label = axis.y.label,
                   legend.title = legend.title,  legend.label = legend.label, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                   legend.text.size = legend.text.size)
  }
  else{
    autoplot.gmwm1(object, background = background, CI = CI, transparence = transparence, bw = bw, 
                  CI.color = CI.color, line.type = line.type, line.color = line.color,
                  point.size = point.size, point.shape = point.shape,
                  title = title, title.size= title.size, 
                  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size , 
                  axis.x.label = axis.x.label,
                  axis.y.label = axis.y.label,
                  legend.title = legend.title,  legend.label = legend.label, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                  legend.text.size = legend.text.size )
  }
}


#' @title Graph Solution of the Generalized Method of Wavelet Moments for Each Process
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM for each latent process.
#' @method autoplot gmwm2
#' @export
#' @keywords internal
#' @param object A \code{GMWM} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @template CommonParams
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM for each latent process.
#' @author JJB, Wenchao
autoplot.gmwm2 = function(object, CI = T, background = 'white', transparence = 0.1, bw = F, 
                          CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                          point.size = NULL,point.shape = NULL,
                          title = NULL, title.size= 15, 
                          axis.label.size = 13, axis.tick.size = 11, 
                          axis.x.label = expression(paste("Scale ", tau)),
                          axis.y.label = expression(paste("Wavelet Variance ", nu)),
                          legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                          legend.text.size = 13,...){
  #require package: grid
  .x=low=high=trans_breaks=trans_format=math_format=value=variable=process=NULL
  
  # Find number of latent processes
  L = length(object$model$desc) + 1
  
  # Get names of latent processes
  # Avoids file system naming issue (e.g. image1, image22, image3)
  nlen = nchar(L)
  nom = sprintf(paste0("z%0",nlen,"d"),1:L)
  
  Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999")
  modulus = (L-1)%/% 9
  remainder = (L-1)%% 9
  process.color = c( rep(Set1, times = modulus), Set1[1:remainder] )
  process.bw_color = gray.colors(L-1, start = 0.2, end = max( max(1, L-2)/(L-1), 0.7))
  
  process.label1 = rep(' ', L-1)
  for(i in 1: (L-1)){
    process.label1[i] = paste0(object$model$desc[i], paste0(rep(' ',times = i), collapse = ''))
  }
  process.label2 = paste0(object$model$desc, collapse = '+')
  process.label = c(process.label1, process.label2)
  
  if(CI == T){
    if(is.null(line.type)){line.type = c('solid','dotted', rep('solid', L) )}
    if(length(line.type)==L+2){line.type = c(line.type[1:2], line.type[2], tail(line.type,L))}
    
    if(is.null(line.color)){line.color = c( "#003C7D", "#999999", process.color, "#F47F24")}
    if(bw){
      line.color = c("#b2b2b2", "#404040", process.bw_color, "#000000")
      CI.color = "grey50"
    }
    if(length(line.color)==L+2){line.color = c(line.color[1:2], line.color[2], tail(line.color,L) )}
    
    if(is.null(point.size)){point.size = c(5,0,rep(0,L-1),5)}
    if(length(point.size)==L+2){point.size = c(point.size[1:2], point.size[2], tail(point.size,L) )}
    
    if(is.null(point.shape)){point.shape = c(20, 46, rep(46,L-1), 1) }
    if(length(point.shape)==L+2){point.shape = c(point.shape[1:2], point.shape[2], tail(point.shape,L) )}
    
    if(is.null(legend.label)){
      #legend.label = c(expression(paste("Empirical WV ", hat(nu))), 
      #                                         expression(paste("CI(", hat(nu)," , 0.95)" )),
      #                                        process.label)
      legend.label = c(bquote("Empirical WV"~hat(nu)), 
                       bquote("CI("*hat(nu)*", "*.(1 - object$alpha)*")" ),
                       process.label) 
    }
    
    df = data.frame(scale = rep(object$scales,L), WV = c(as.vector(object$decomp.theo), object$theo), 
                    process = rep(nom, each = length(object$scales)))
    
    WV = data.frame(emp = object$wv.empir, low = object$ci.low, high = object$ci.high, scale = object$scales)
    melt.wv = melt(WV, id.vars = 'scale')
    
    breaks = c('emp','low',nom)
    legend.fill = c(NA, CI.color, rep(NA, L) )
    legend.linetype = c(line.type[1], 'blank', line.type[4:length(line.type)])
    legend.pointshape = c(point.shape[1],NA, point.shape[4:length(point.shape)])
    ##legend.pointsize = c(point.size[1:2],0) DON'T NEED THIS LINE
    
  }else{
    if(is.null(line.type)){line.type = c('solid', rep('solid', L))}
    
    if(is.null(line.color)){line.color = c( "#003C7D", process.color, "#F47F24")}
    if(bw){
      #line.color = c("#000000", "#b2b2b2", "#404040")
      line.color = c("#b2b2b2", process.bw_color, "#000000" )
    }
    
    if(is.null(point.size)){point.size = c(5, rep(0,L-1), 5)}
    if(is.null(point.shape)){point.shape =c(20, rep(46,L-1), 1 ) }
    
    if(is.null(legend.label)){legend.label = c(expression(paste("Empirical WV ", hat(nu))), 
                                               process.label)}
    
    df = data.frame(scale = rep(object$scales,L), WV = c(as.vector(object$decomp.theo), object$theo), 
                    process = rep(nom, each = length(object$scales)))
    WV = data.frame(emp = object$wv.empir, scale = object$scales)
    melt.wv = melt(WV, id.vars = 'scale')
    
    breaks = c('emp', nom)
    #legend.fill = c(rep(NA, L+1))
    #legend.linetype = c(line.type[1:(L+1)],'blank')
    #legend.pointshape = c(point.shape[1:(L+1)],NA)
  }
  
  p = ggplot() + 
    geom_line( data = melt.wv, mapping = aes(x = scale, y = value, color = variable, linetype = variable) )+
    geom_point(data = melt.wv, mapping = aes(x = scale, y = value, color = variable, size = variable, shape = variable))
  
  p = p + geom_line(data = df, mapping = aes(x = scale, y = WV,color = process, linetype = process)) + 
    geom_point(data = df, mapping = aes(x = scale, y = WV, color = process, size = process, shape = process)) + 
    
    scale_linetype_manual(name = legend.title, values = c(line.type),breaks = breaks, labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label)+
    scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label)
  
  if(CI){
    p = p +
      geom_ribbon(data = WV, mapping = aes(x = scale, y = NULL,ymin =low, ymax = high ), fill = CI.color, alpha = transparence, show_guide = T) +
      #scale_fill_manual(name = legend.title, values = c(color.CI,'red'), breaks = breaks, labels = legend.label) +
      guides(colour = guide_legend(override.aes = list(fill = legend.fill, linetype = legend.linetype, shape = legend.pointshape)))
  }
  
  p = p +  
    scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x)) ) +
    
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    coord_cartesian(ylim = c(0.8* min( c(object$ci.low, object$wv.empir) ), 1.05*max( c(object$wv.empir, object$ci.high))) )
  
  if( background == 'white' || bw){
    p = p + theme_bw() 
  }
  
  p = p +
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      legend.text.align = 0 )
    
  
  if (is.null(title)){
    if(object$robust){
      p = p + ggtitle("Haar Wavelet Variance Robust Representation")
    }
    else{
      p = p + ggtitle("Haar Wavelet Variance Classical Representation")
    }
  }
  
  p
}


#' @title Graph Solution of the Generalized Method of Wavelet Moments Non-individually
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method autoplot gmwm1
#' @export
#' @keywords internal
#' @param object A \code{GMWM} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @template CommonParams
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB, Wenchao
autoplot.gmwm1 = function(object, CI = T, background = 'white', transparence = 0.1, bw = F, 
                         CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                         point.size = NULL, point.shape = NULL,
                         title = NULL, title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)),
                         legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  #require pakage: scales, grid
  low=high=emp=theo=trans_breaks=trans_format=math_format=.x=value=variable=NULL

  temp = data.frame( emp = object$wv.empir,
                     low = object$ci.low,
                     high = object$ci.high,
                     theo = object$theo,
                     scale = object$scales)
  
  if(CI == T){
    if(is.null(line.type)){line.type = c('solid', 'dotted', 'solid')}
    if(length(line.type)==3){line.type = c(line.type[1:2], line.type[2:3])}
    
    if(is.null(line.color)){line.color = c("#003C7D", "#999999" , "#F47F24")}
    if(bw){
      line.color = c("#b2b2b2", "#404040", "#000000")
      color.CI = "grey50"
    }
    if(length(line.color)==3){line.color = c(line.color[1:2],line.color[2:3])}
    
    if(is.null(point.size)){point.size = c(5, 0, 5)}
    if(length(point.size)==3){point.size = c(point.size[1:2],point.size[2:3])}
    
    if(is.null(point.shape)){point.shape = c(20, 46, 1) }
    if(length(point.shape)==3){point.shape = c(point.shape[1:2],point.shape[2:3])}
    
    if(is.null(legend.label)){
      #legend.label = c(expression(paste("Empirical WV ", hat(nu))), 
      #                                         expression(paste("CI(", hat(nu)," , 0.95)" )),
      #                                         expression(paste("Implied WV ", nu,"(",hat(theta),")")) )
      
      legend.label = c(bquote("Empirical WV"~hat(nu)), 
                       bquote("CI("*hat(nu)*", "*.(1 - object$alpha)*")" ),
                       bquote("Implied WV"~nu*"("*hat(theta)*")")) 
    }
    
    WV = melt(temp, id.vars = 'scale')
    breaks = c('emp','low','theo')
    legend.fill = c(NA,CI.color,NA )
    legend.linetype = c(line.type[1],'blank', tail(line.type,1) )
    legend.pointshape = c(point.shape[1], NA, tail(point.shape,1) )
    #legend.pointsize = c(point.size[1:2],0)
  }else{
    if(is.null(line.type)){line.type = c('solid','solid')}
    if(is.null(line.color)){line.color = c("#003C7D", "#F47F24")}
    if(bw){
      line.color = c("#b2b2b2", "#000000")}
    if(is.null(point.size)){point.size = c(5,5)}
    if(is.null(point.shape)){point.shape = c(20,1) }
    if(is.null(legend.label)){legend.label = c(expression(paste("Empirical WV ", hat(nu))),
                                               expression(paste("Implied WV ", nu,"(",hat(theta),")"))    )}
    
    WV = melt(temp, id.vars = 'scale', measure.vars = c('emp', 'theo'))
    breaks = c('emp','theo')
    #legend.color = c(NA,NA)
  }
  
  p = ggplot(data = WV, mapping = aes(x = scale)) + geom_line(aes(y = value, color = variable, linetype = variable)) +
    geom_point(aes(y = value, shape = variable, size = variable, color = variable)) + 
    scale_linetype_manual(name = legend.title, values = c(line.type), breaks = breaks, labels = legend.label) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks,labels = legend.label)+
    
    scale_size_manual(name = legend.title, values = c(point.size),breaks = breaks,labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label) 
  
#   obj2 = data.frame(y = WV$value[WV$variable == 'theo'], x = WV$scale[WV$variable == 'theo'])
#   p = p + geom_point(data = obj2, mapping = aes(x = x, y = y), size = 5, shape = 1, colour = line.color[2]) 
#   p = p + geom_point(data = obj2, mapping = aes(x = x, y = y), size = 4.5, shape = 1, colour = line.color [2])
  
  if(CI){
    p = p +
      geom_ribbon(data = temp, mapping = aes(ymin = low, ymax = high),fill = CI.color, show_guide = T,alpha = transparence) +
      
      #scale_fill_manual(name = legend.title, values = c(color.CI,'red'), breaks = breaks, labels = legend.label) +
      guides(colour = guide_legend(override.aes = list(fill = legend.fill, linetype = legend.linetype, shape = legend.pointshape)))
  }
  
  #decide where to place the legend
  legendPlace = placeLegend(temp$emp[1], temp$low[ length(temp$low) ], temp$high[ length(temp$high)])    
  p = p + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  if( background == 'white' || bw){
    p = p + theme_bw() 
  }
  
  p = p +
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      legend.text.align = 0 )
  # if(!bw){p = p + theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))}
  if (is.null(title)){
    if(object$robust){
      p = p + ggtitle("Haar Wavelet Variance Robust Representation")
    }
    else{
      p = p + ggtitle("Haar Wavelet Variance Classical Representation")
    }
  }
  p
}


#' @title Graphically Compare GMWM Model Fit
#' @description Creates GMWM model fits of different models within the same panel.
#' @param ... Several \code{gmwm} objects
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE).
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param auto.label.wvar A \code{boolean} that indicates whether legend label should indicate the gmwm objects are robust or classical
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of the confidence interval (e.g. 'black', 'red', '#003C7D', etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param legend.title A \code{string} that indicates the title of legend.
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend.
#' @param nrow An \code{integer} that indicates how many rows the graph should be arranged in.
#' @param plot.emp.wv A \code{boolean} that indicates whether Emp. WV should be plotted or not (Used in \code{compare.models}).
#' @return A ggplot2 panel containing the graph of gmwm objects.
#' @author JJB, Wenchao
#' @details 
#' If only one object is supplied, this function is actually calling \code{plot.gmwm}. When the parameters \code{line.color}, 
#' \code{CI.color}, \code{line.type},  \code{point.size}, \code{point.shape} and \code{legend.label} are modified, please follow the rules of \code{\link{plot.gmwm}}.
#' 
#' When \code{CI = T}, for \code{CI.color}, specify the color for each object. 
#' For \code{line.color}, \code{line.type}, \code{point.size}, \code{point.shape},
#' specify the value of lower bound, upper bound, empirical wavelet variance (WV), implied WV respectively for each object.
#' 
#' When \code{CI = F}, you don't need \code{CI.color} this time. For \code{line.color}, \code{line.type}, \code{point.size}, \code{point.shape},
#' only specify the value of empirical WV and implied WV respectively.
#' 
#' Check the examples for help.
#' 
#' @examples
#' \dontrun{# AR
#' set.seed(8836)
#' n = 200
#' x = gen.gts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.gmwm(GMWM1, GMWM2, split = FALSE)
#' compare.gmwm(GMWM1, GMWM2, point.size = rep(c(1,1,4,4),2), CI.color = c('black','grey'))
#' compare.gmwm(GMWM1, GMWM2, CI = F, point.size = rep(c(6,6),2))
#' }
compare.gmwm = function(..., background = 'white', split = TRUE, CI = TRUE, auto.label.wvar = F, transparence = 0.1, line.color = NULL, 
                        CI.color = NULL, line.type = NULL,  point.size = NULL, point.shape = NULL,
                        title = "Comparison of Implied Wavelet Variance", title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu)),
                        facet.label.size = 13,
                        legend.label = NULL,
                        legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                        legend.text.size = 13, nrow = 1, plot.emp.wv = T ){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...)
  numObj = length(obj_list)
  object.names = as.character(substitute(...()))
  if(any( count_models(object.names) >1 )){
    stop('Duplicated objects are detected.')
  }
  
  if (numObj == 0){
    stop('At least one gmwm object should be given')
  }else if(numObj >2){
    stop('The function can only deal with at most 2 objects currently. This constraint may be removed in next version.')
  }else if (numObj == 1){
    ## just plot
    if(is.null(CI.color)){CI.color = "#003C7D"}
    warning('One object is supplied. You are actually calling plot.gmwm().')
    
    plot(..., process.decomp = FALSE, background = background,
         CI = CI, transparence = transparence, bw = F, CI.color = CI.color,
         line.type = line.type, line.color = line.color, point.size = point.size,
         point.shape = point.shape, title = title, title.size = title.size, axis.label.size = axis.label.size,
         axis.tick.size = axis.tick.size, axis.x.label = axis.x.label,
         axis.y.label = axis.y.label,
         legend.title = legend.title, legend.label = legend.label, legend.key.size = legend.key.size,
         legend.title.size = legend.title.size, legend.text.size = legend.text.size)
  }
  else  {
    
    #check parameter
    
    if(CI){
      params = c('line.color', 'CI.color', 'line.type', 'point.size', 'point.shape', 'legend.label')
      requireLength = c(4*numObj, numObj, 4*numObj, 4*numObj, 4*numObj, 2*numObj)
      default = list(NULL, NULL,  NULL, NULL, NULL, NULL)
      nullIsFine = c(rep(T,6))
    }else{
      params = c('line.color', 'line.type', 'point.size', 'point.shape', 'legend.label')
      requireLength = c(2*numObj, 2*numObj, 2*numObj, 2*numObj, 2*numObj)
      default = list(NULL, NULL,  NULL, NULL, NULL)
      nullIsFine = c(rep(T,5))
    }
    
    for (i in 1:length(params)){
      one_param = params[i]
      if( length(get(one_param))!=requireLength[i]){
        isNull = is.null(get(one_param))
        if(isNull && nullIsFine[i]){}else{
          warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                        'is supplied.','Default setting is used.'))
        }
        assign(one_param, default[[i]])
      }
    }
    
    if(CI){
      #supply order: low, high, emp.wv, imp.wv --- rep it for numObj times
      if(is.null(point.size)){point.size = rep(c(0,0,5,5), numObj)}
      if(is.null(point.shape)){point.shape = rep(c(20, 20, 20, 1), numObj)}
      if(is.null(line.type)){line.type = rep(c('dotted','dotted', 'solid','solid'), numObj)}
      
      if(is.null(line.color)){
        if (numObj == 2){
          wv.palette = c("#003C7D","#F47F24")
        }else{
          # 'Dark2'
          Dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
          modulus = numObj%/% 8
          remainder = numObj%% 8
          wv.palette = c( rep(Dark2, times = modulus), Dark2[1:remainder] )
        }
        theo.palette = ggColor(numObj)
        
        line.color = rep(NA, numObj* 4) #initialize
        for (i in 1:numObj){
          for (j in 1:4){
            if(j==1||j==2||j==3){
              line.color [4*i-4 + j] = wv.palette[i]
            }else{
              line.color[4*i] = theo.palette[i]
            }
          }
        }
        
        # if attributes are same, use same color
        for(i in 1:(numObj-1)){
          for(j in (i+1):numObj ){
            same.scales = length(obj_list[[i]]$scales) == length(obj_list[[j]]$scales)
            if(same.scales){
              if(all(obj_list[[i]]$scales == obj_list[[j]]$scales)){
                same.scales = T
              }else{same.scales = F}
            }
            
            if(same.scales){
              same.wv = obj_list[[i]]$expect.diff == obj_list[[j]]$expect.diff
              #same.wv = all(obj_list[[i]]$wv.empir == obj_list[[j]]$wv.empir)
              same.theo = length(obj_list[[i]]$theo) == length(obj_list[[j]]$theo)
              if(same.theo){
                if(all(obj_list[[i]]$theo == obj_list[[j]]$theo)){
                  same.theo = T
                }else{
                  same.theo = F
                }
              }
              
              same.alpha = obj_list[[i]]$alpha == obj_list[[j]]$alpha
              if(same.wv){line.color[4*(j-1)+3] = line.color[4*(i-1)+3]}
              if(same.theo){line.color[4*(j-1)+4] = line.color[4*(i-1)+4]}
              if(same.wv && same.alpha){
                line.color[4*(j-1)+1] = line.color[4*(i-1)+1]
                line.color[4*(j-1)+2] = line.color[4*(i-1)+2]
              }
            }
          }
        }#end: attributes check
        
      }
      
      if(is.null(CI.color)){
        CI.color = rep(NA, numObj)
        for(i in 1:numObj){
          CI.color[i] = line.color[4*(i-1) +1]
        }
        CI.color = alpha(CI.color, transparence)
      }
      
      #for legend.label
      #wv_string = expression(paste("Empirical WV ", hat(nu),' and ','CI') )
      wv_string = bquote("Empirical WV"~hat(nu)~'and'~'CI')
    }else{
      if(is.null(point.size)){point.size = rep(c(5,5), numObj)}
      if(is.null(point.shape)){point.shape = rep(c(20, 1), numObj)}
      if(is.null(line.type)){line.type = rep(c('solid','solid'), numObj)}
      
      if(is.null(line.color)){
        
        #wv.palette = ggColor(numObj)
        #if (numObj == 2){
        #  wv.palette = c("#003C7D","#F47F24")}
        #theo.palette = alpha(wv.palette, 0.7)
        
        if (numObj == 2){
          wv.palette = c("#003C7D","#F47F24")
        }else{
          # 'Dark2'
          Dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
          modulus = numObj%/% 8
          remainder = numObj%% 8
          wv.palette = c( rep(Dark2, times = modulus), Dark2[1:remainder] )
        }
        theo.palette = ggColor(numObj)
        
        
        line.color = rep(NA, numObj* 2) #initialize
        for (i in 1:numObj){
          for (j in 1:2){
            if(j==1){
              line.color [2*(i-1) + j] = wv.palette[i]
            }else{
              line.color[2*i] = theo.palette[i]
            }
          }
        }
        
        # if attributes are same, use same color
        for(i in 1:(numObj-1)){
          for(j in (i+1):numObj ){
            same.scales = length(obj_list[[i]]$scales) == length(obj_list[[j]]$scales)
            if(same.scales){
              if(all(obj_list[[i]]$scales == obj_list[[j]]$scales)){
                same.scales = T
              }else{same.scales = F}
            }
            
            if(same.scales){
              same.wv = obj_list[[i]]$expect.diff == obj_list[[j]]$expect.diff
              #same.wv = all(obj_list[[i]]$wv.empir == obj_list[[j]]$wv.empir)
              same.theo = length(obj_list[[i]]$theo) == length(obj_list[[j]]$theo)
              if(same.theo){
                if(all(obj_list[[i]]$theo == obj_list[[j]]$theo)){
                  same.theo = T
                }else{
                  same.theo = F
                }
              }
              
              if(same.wv){line.color[2*(j-1)+1] = line.color[2*(i-1)+1]}
              if(same.theo){line.color[2*(j-1)+2] = line.color[2*(i-1)+2]}
            }
          }
        }#end: attributes check
        
      }
      
      #for legend.label
      #wv_string = expression(paste("Empirical WV ", hat(nu)))
      wv_string = bquote("Empirical WV" ~hat(nu))
    }#end CI
    
    #legend.label
    #theo_string = expression(paste("Implied WV ", nu,"(",hat(theta),")"))
    theo_string = bquote("Implied WV"~nu*"("*hat(theta)*")")
    
    if(is.null(legend.label)){
      legend.label = c()
      legend.label_raw = c()
      for (i in 1:(2* numObj) ){
        legend.label_raw[i] = paste( object.names[(i-1)%/%2 + 1] )
        if(auto.label.wvar){
          legend.label_raw[i] = paste(legend.label_raw[i], if(obj_list[[(i-1)%/%2 + 1]]$robust) '(Robust)' else '(Classical)')}
        if(i%%2 == 1){
          legend.label[i] = as.expression(bquote(.(legend.label_raw[i])~ .(wv_string) ) )
        }else{
          legend.label[i] = as.expression(bquote( .(legend.label_raw[i])~ .(theo_string)) )
        }
      }
    }
    
    #breaks
    breaks = rep(NA, 2*numObj)
    for(i in 1:(2* numObj)){
      if(i%%2 == 1){
        breaks[i] = paste( object.names[(i-1)%/%2 + 1], 'WV')
      }else{
        breaks[i] = paste( object.names[(i-1)%/%2 + 1], 'z_theo')
      }
    }
    
    #levels
    if(plot.emp.wv){#if plot WV
      levels = rep(NA, 4*numObj)
      for(i in 1:(4*numObj)){
        if(i%%4 == 1){
          levels[i] = paste(object.names[(i-1)%/%4 + 1], 'low')
        }else if(i%%4 == 2){
          levels[i] = paste(object.names[(i-1)%/%4 + 1], 'high')
        }else if(i%%4 == 3){
          levels[i] = paste(object.names[(i-1)%/%4 + 1], 'WV')
        }else{
          levels[i] = paste(object.names[(i-1)%/%4 + 1], 'z_theo')
        }
      }
    }else{#if not plot WV, only used in compare.models. At the same time, CI = F
      levels = rep(NA, numObj)
      for(i in 1:numObj){
        levels[i] = paste(object.names[i], 'z_theo')
      }
    }
    
    
    total.len = 0
    each.len = numeric(numObj)
    for (i in 1:numObj){
      each.len[i] = length(obj_list[[i]]$wv.empir)
      total.len = total.len + each.len[i]
    }
    
    if(CI){
      #Initialize empty data frame with right number of rows
      obj = data.frame(WV = numeric(total.len),
                       scales = numeric(total.len),
                       low = numeric(total.len),
                       high = numeric(total.len),
                       z_theo = numeric(total.len),
                       dataset = 'XYZ', stringsAsFactors=FALSE)
      
      #put data into data frame
      t = 1
      for (i in 1:numObj){
        d = each.len[i]
        obj[t:(t+d-1),] = data.frame(WV = obj_list[[i]]$wv.empir,
                                     scales = obj_list[[i]]$scales,
                                     low = obj_list[[i]]$ci.low,
                                     high = obj_list[[i]]$ci.high,
                                     z_theo = obj_list[[i]]$theo,
                                     dataset = object.names[i], stringsAsFactors=FALSE)
        t = t +d
      }
    }else{
      #Initialize empty data frame with right number of rows
      obj = data.frame(WV = numeric(total.len),
                       scales = numeric(total.len),
                       z_theo = numeric(total.len),
                       dataset = 'XYZ', stringsAsFactors=FALSE)
      
      #put data into data frame
      t = 1
      for (i in 1:numObj){
        d = each.len[i]
        obj[t:(t+d-1),] = data.frame(WV = obj_list[[i]]$wv.empir,
                                     scales = obj_list[[i]]$scales,
                                     z_theo = obj_list[[i]]$theo,
                                     dataset = object.names[i], stringsAsFactors=FALSE)
        t = t +d
      }
    }
    
    melt.obj = melt(obj, id.vars = c('scales', 'dataset'))
    #if (numObj == 2 ){
    autoplot.gmwmComp(melt.obj, breaks = breaks, split = split, CI = CI, background = background, transparence = transparence, line.color =line.color, 
                      CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
                      title = title, title.size= title.size, 
                      axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                      axis.x.label = axis.x.label,
                      axis.y.label = axis.y.label,
                      facet.label.size = facet.label.size,
                      legend.label = legend.label,
                      legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                      legend.text.size = legend.text.size,
                      nrow = nrow, plot.emp.wv = plot.emp.wv, object.names = object.names, levels = levels)
    #     }
    #     else{
    #       #will be modified in next version
    #       autoplot.gmwmComp(melt.obj, breaks = breaks, split = split, CI = CI, background = background, transparence = transparence, line.color = line.color, 
    #                         CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
    #                         title = title, title.size= title.size, 
    #                         axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
    #                         axis.x.label = axis.x.label,
    #                         axis.y.label = axis.y.label,
    #                         facet.label.size = facet.label.size,
    #                         legend.label = legend.label,
    #                         legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
    #                         legend.text.size = legend.text.size,
    #                         nrow = nrow, plot.emp.wv = plot.emp.wv)
    #     }
    
  }
  
}

#' @title Compare GMWM Model Fits with ggplot2 (Internal)
#' @description Creates a single graph that contains several GMWM models plotted against each other.
#' @method autoplot gmwmComp
#' @export
#' @keywords internal
#' @param object A \code{data.frame} containing both sets of GMWM object data.
#' @param breaks A \code{vector} used to determine the legend label.
#' @param levels A \code{vector} of \code{string} that indicates each level in the dataset.
#' @param object.names A \code{vector} of \code{string} that indicates name of each object.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE).
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of the confidence interval (e.g. 'black', 'red', '#003C7D', etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param legend.title A \code{string} that indicates the title of legend.
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend.
#' @param nrow An \code{integer} that indicates how many rows the graph should be arranged in.
#' @param plot.emp.wv A \code{boolean} that indicates whether Emp. WV should be plotted or not (Used in \code{compare.models}).
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing one graph with several GMWM models plotted against each other.
#' @note User doesn't need to know this function.
#' @author JJB, Wenchao
autoplot.gmwmComp = function(object, breaks, levels, object.names, split = TRUE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             CI.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                             title = "Comparison of Implied Wavelet Variance", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)),
                             facet.label.size = 13,
                             legend.label = NULL,
                             legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                             legend.text.size = 13, nrow = 1, plot.emp.wv = T, ...){
  scales=low=high=WV=emp=theo=trans_breaks=trans_format=math_format=.x=dataset=value=variable=NULL
  
  if(CI){object.CI = object[object$variable =='low'|object$variable=='high', ]}
  if(!plot.emp.wv && !CI){
    object = object[object$variable != 'WV', ]
    line.color = line.color[seq(from = 2, to = length(line.color), by = 2)]
    line.type = line.type[seq(from = 2, to = length(line.type), by = 2)]
    point.size = point.size[seq(from = 2, to = length(point.size), by = 2)]
    point.shape = point.shape[seq(from = 2, to = length(point.shape), by = 2)]
    legend.label = legend.label[seq(from = 2, to = length(legend.label), by = 2)]
    breaks = breaks[seq(from = 2, to = length(breaks), by = 2)]
  }
  object$variable = paste(object$dataset, object$variable)
  object$dataset = factor(object$dataset, levels = object.names)
  object$variable = factor(object$variable, levels = levels)
  
  p = ggplot() + 
    geom_line( data = object, mapping = aes(x = scales, y = value, color = variable, linetype = variable)) + 
    geom_point(data = object, mapping = aes(x = scales, y = value, color = variable, size = variable, shape = variable)) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  p = p + scale_size_manual(name = legend.title, values = point.size, breaks = breaks,labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = point.shape, breaks = breaks,labels = legend.label) +
    scale_linetype_manual(name = legend.title, values = line.type, breaks = breaks,labels = legend.label) +
    scale_color_manual(name = legend.title, values = line.color, breaks = breaks,labels = legend.label)
  
  if(CI){
    object.CI = dcast(object.CI, scales+ dataset~variable)
    object.CI$dataset = factor(object.CI$dataset, levels = object.names )
    
    p = p + 
      geom_ribbon(data = object.CI, mapping = aes(x = scales, ymin = low, ymax = high, fill = dataset), alpha = transparence, show_guide = T)
    
    legend.fill = rep(NA, 2*length(CI.color))
    for(i in 1:length(legend.fill)){
      if(i%%2 == 1){
        legend.fill[i] = CI.color[(i-1)%/%2+1]
      }else{
        legend.fill[i] = NA
      }
    }
    
    #get CI.color from above function
    p = p + 
      scale_fill_manual(name = legend.title, values = CI.color, breaks = breaks,labels = legend.label) +
      guides(colour = guide_legend(override.aes = list(fill = legend.fill)) )
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  if (split){
    p = p + facet_wrap(~dataset,nrow = nrow) 
    #+ theme(legend.position="none")
  }
  
  p = p +  xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      strip.text = element_text(size = facet.label.size),
      legend.text.align = 0) 
  
  p
}




#' @title Graphically Compare GMWM Models Constructed by the Same Data
#' @description Creates a table of graphs to compare GMWM model fits.
#' @param ... Several \code{gmwm} objects, and they must be constrcuted by the same data.
#' @param display.model A \code{boolean} indicating whether the model should be displayed in the facet label.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of confidence interval.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. 'black', 'red', '#003C7D', etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @author Wenchao
#' @details 
#' This function only works for \code{gmwm} objects which are constrcuted by same data, and 
#' all \code{gmwm} objects must be constructed by classical method, or by robust methods
#' with the same efficiency. That's because this function assumes each \code{gmwm} object has the same empirical wavelet variance (WV).
#' This function will check whether this requirement is satisfied before plotting the graph.
#' 
#' \code{line.type}, \code{line.color}, \code{point.size}, \code{point.size} must be a \code{vector}. You need to follow this order:
#' "WV, bounds of CI, model1 implied WV, model2 implied WV, ...."
#' 
#' Please check examples for help.
#' 
#' If you meet the error "polygon edge not found", it is complaining that you don't have enough space to
#' plot the graph. Adjust the plot window.
#' @examples
#' \dontrun{
#' library(GMWM)
#' data(imu)
#' model1 = gmwm.sensor(3*AR1(),imu[,2])
#' model2 = gmwm.sensor(2*AR1() + RW(),imu[,2])
#' compare.models(model1, model2)
#' compare.models(model1, model2, display.model = F, point.size = c(4, 0, 4, 4))
#' compare.models(model1, model2, transparence = 0.2, line.color = c('black', 'grey', 'blue', 'red'))
#' }
compare.models = function(..., display.model = T, background = 'white', transparence = 0.1, CI.color = "#003C7D",
                          line.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                          title = "Comparison of GMWM Models", title.size= 18, 
                          axis.label.size = 16, axis.tick.size = 11, 
                          facet.label.size = 13,
                          axis.x.label = expression(paste("Scale ", tau)),
                          axis.y.label = expression(paste("Wavelet Variance ", nu))){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...)
  numObj = length(obj_list)
  
  for(i in 1:numObj){
    if( !is(obj_list[[i]], 'gmwm') ){
      stop('Model you supplied was not gmwm object.')
    }
  }
  
  if(is.null(CI.color)){
    CI.color = "#003C7D"
  }
  
  if(numObj<1){
    stop('At least one model must be supplied.')
  }
  object.names = as.character(substitute(...()))
  
  #check parameter
  params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape')
  requireLength = c(2+numObj, 2+numObj, 1, 2+numObj, 2+numObj)
  default = list(NULL, NULL,  "#003C7D", NULL, NULL)
  nullIsFine = c(rep(T,5))
  for (i in 1:length(params)){
    one_param = params[i]
    if( length(get(one_param))!=requireLength[i]){
      isNull = is.null(get(one_param))
      if(isNull && nullIsFine[i]){}else{
        warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                      'is supplied.','Default setting is used.'))
      }
      assign(one_param, default[[i]])
    }
  }
  
  if(numObj == 1){
    #plot.gmwm will do the parameter checking
    #other parameters are not listed here. But they cannot be passed to plot.gmww by '...'
    warning('One object is supplied. You are actually calling plot.gmwm().')
    plot.gmwm(x = obj_list[[1]], process.decomp = F, background = background, transparence = transparence, CI.color = CI.color, CI = T, bw = F, line.type = line.type,
              line.color = line.color, point.size = point.size, point.shape = point.shape, title = title,
              title.size = title.size, axis.label.size = axis.label.size, axis.tick.size = axis.tick.size,
              axis.x.label = axis.x.label, axis.y.label = axis.y.label)
    
  }else{
    
    for (i in 1:(numObj-1) ){
      #1. check expect.diff
      if ( obj_list[[i]]$expect.diff != obj_list[[i+1]]$expect.diff ){
        stop('This function can only operate on models constrcuted by the same data.') 
      }
      
      if(!obj_list[[1]]$robust){
        #2.1 classical case, make sure all other objects are classical
        if(obj_list[[i+1]]$robust != F){
          stop('Make sure your models are all classical, or they are all robust and with same efficiency.')
        }
      }else{
        #2.2 robust case, make sure all other objects are robust and have same eff
        if(obj_list[[i+1]]$robust != T || (obj_list[[i]]$eff != obj_list[[i+1]]$eff)){
          stop('Make sure your models are all classical, or they are all robust and with same efficiency.')
        }
      }
    }
    
    #create one empty list
    plot.list = vector('list', numObj * numObj)
    
    #supply order: wv, high/low, model1 imp.wv, model2 imp.wv, model3 imp.wv ....
    if(is.null(line.color)){
      wv.color = "#003C7D" #color for emp.wv
      theo.color = ggColor(numObj) #color for imp.wv
      high_low.color = wv.color
    }else{
      wv.color = line.color[1] #color for emp.wv
      theo.color = line.color[3:length(line.color)] #color for imp.wv
      high_low.color = line.color[2]
    }
    #line.color1 is supplied to plot.gmwm
    line.color1 = c(wv.color, high_low.color, theo.color)
    
    #line.color2 is supplied to compare.gmwm when CI = T
    line.color2 = c(rep(high_low.color, 2), wv.color, theo.color)
    
    #line.color3 is supplied to compare.gmwm when CI = F
    line.color3 = c(wv.color, theo.color)
    
    
    if(is.null(line.type)){
      wv.line.type = "solid" # for emp.wv
      theo.line.type = rep('solid', numObj) # for imp.wv
      high_low.line.type = 'dotted'
      
    }else{
      wv.line.type = line.type[1]
      theo.line.type = line.type[3:length(line.type)]
      high_low.line.type = line.type[2]
    }
    #line.type1 is supplied to plot.gmwm
    line.type1 = c(wv.line.type, high_low.line.type, theo.line.type)
    
    #line.type2 is supplied to compare.gmwm when CI = T
    line.type2 = c(rep(high_low.line.type, 2), wv.line.type, theo.line.type)
    
    #line.type3 is supplied to compare.gmwm when CI = F
    line.type3 = c(wv.line.type, theo.line.type)
    
    
    if(is.null(point.size)){
      if(numObj > 10){temp.point.size = 1
      }else{
        temp.point.size = 4.7-0.4*(numObj-1)}
      
      wv.point.size = temp.point.size # for emp.wv
      theo.point.size = rep(temp.point.size, numObj) +0.1 # for imp.wv
      high_low.point.size = 0
    }else{
      wv.point.size = point.size[1]
      theo.point.size = point.size[3:length(point.size)]
      high_low.point.size = point.size[2]
    }
    #point.size1 is supplied to plot.gmwm
    point.size1 = c(wv.point.size, high_low.point.size, theo.point.size)
    
    #point.size2 is supplied to compare.gmwm when CI = T
    point.size2 = c(rep(high_low.point.size, 2), wv.point.size, theo.point.size)
    
    #point.size3 is supplied to compare.gmwm when CI = F
    point.size3 = c(wv.point.size, theo.point.size)
    
    
    if(is.null(point.shape)){
      wv.point.shape = 20 # for emp.wv
      theo.point.shape = rep(1, numObj) # for imp.wv
      high_low.point.shape = 20
      
    }else{
      wv.point.shape = point.shape[1]
      theo.point.shape = point.shape[3:length(point.shape)]
      high_low.point.shape = point.shape[2]
    }
    #point.shape1 is supplied to plot.gmwm
    point.shape1 = c(wv.point.shape, high_low.point.shape, theo.point.shape)
    
    #point.shape2 is supplied to compare.gmwm when CI = T
    point.shape2 = c(rep(high_low.point.shape, 2), wv.point.shape, theo.point.shape)
    
    #point.shape3 is supplied to compare.gmwm when CI = F
    point.shape3 = c(wv.point.shape, theo.point.shape)
    
    
    #CI.color2 is supplied to compare.gmwm when CI = T
    CI.color2 = rep(CI.color, 2)
    
    
    #common setting
    commonSet = theme(legend.position = 'none') #remove legend
    #axis.title = element_blank(), #remove the space for axes
    #plot.margin = unit(c(0,0,0,0), "mm") )#remove margin around the graph
    
    #keep x
    keep.x = theme(axis.text.x=element_text(size=axis.tick.size),
                   axis.title.x=element_text(size=facet.label.size))
    #axis.ticks.x=element_line(color = "#0000FF00")
    #remove x
    remove.x = theme(axis.text.x=element_text(size=axis.tick.size, color="#0000FF00"),
                     axis.title.x=element_text(size=facet.label.size, color = "#0000FF00"),
                     axis.ticks.x=element_line(color = "#0000FF00") )
    #keep y
    keep.y = theme(axis.text.y=element_text(size=axis.tick.size),
                   axis.title.y=element_text(size=facet.label.size))
    #axis.ticks.y=element_line(color = "#0000FF00") )
    #remove y
    remove.y = theme(axis.text.y=element_text(size=axis.tick.size, color="#0000FF00"),
                     axis.title.y=element_text(size=facet.label.size, color="#0000FF00"),
                     axis.ticks.y=element_line(color = "#0000FF00"))
    #keep title
    keep.title = theme(plot.title = element_text(size=facet.label.size))
    #remove title
    remove.title = theme(plot.title = element_text(size=facet.label.size, color="#0000FF00") )
    
    n = numObj
    for (i in 1:n){ #i: row index
      for (j in 1:n){ #j: col index
        
        #select right title and axis.y.label
        if(display.model){new.title = getModel.gmwm(obj_list[[j]])} else{new.title = object.names[j]}
        if(display.model){new.axis.y.label = getModel.gmwm(obj_list[[i]])} else{new.axis.y.label = object.names[i]}
        new.axis.x.label = ''
        
        if(i == j){
          plot.list[[(i-1)*n+j]] = plot.gmwm( obj_list[[i]], process.decomp = F, background = background,
                                              CI = T, transparence = transparence, bw = F, CI.color = CI.color, 
                                              line.type = c(line.type1[1:2], line.type1[i+2]), line.color = c(line.color1[1:2], line.color1[i+2]), 
                                              point.size = c(point.size1[1:2], point.size1[i+2]), 
                                              point.shape = c(point.shape1[1:2], point.shape1[i+2]), 
                                              title = new.title, title.size = facet.label.size, axis.label.size = facet.label.size,
                                              axis.tick.size = axis.tick.size, axis.x.label = new.axis.x.label, axis.y.label = new.axis.y.label ) + 
            commonSet
          
        }else if(i>j){
          
          #reconstruct line.color
          line.color = c(line.color2[1:3], line.color2[3+i], line.color2[1:3], line.color2[3+j])
          #reconstruct line.type
          line.type = c(line.type2[1:3], line.type2[3+i], line.type2[1:3], line.type2[3+j])
          #reconstruct point.size
          point.size = c(point.size2[1:3], point.size2[3+i], point.size2[1:3], point.size2[3+j])
          #reconstruct point.shape
          point.shape = c(point.shape2[1:3], point.shape2[3+i], point.shape2[1:3], point.shape2[3+j])
          
          plot.list[[(i-1)*n+j]] = compare.gmwm(obj_list[[i]], obj_list[[j]], split = F, CI = T, auto.label.wvar = F,
                                                transparence = transparence/2, line.type = line.type, point.size = point.size,
                                                point.shape = point.shape, CI.color = CI.color2,
                                                line.color = line.color, 
                                                title = new.title, title.size = facet.label.size, 
                                                axis.label.size = facet.label.size, axis.tick.size = axis.tick.size,
                                                axis.x.label = new.axis.x.label, axis.y.label = new.axis.y.label, plot.emp.wv = T ) +
            commonSet
          
        }else{
          #reconstruct line.color
          line.color = c(line.color3[1], line.color3[1+i], line.color3[1], line.color3[1+j])
          #reconstruct line.type
          line.type = c(line.type3[1], line.type3[1+i], line.type3[1], line.type3[1+j])
          #reconstruct point.size
          point.size = c(point.size3[1], point.size3[1+i], point.size3[1], point.size3[1+j])
          #reconstruct point.shape
          point.shape = c(point.shape3[1], point.shape3[1+i], point.shape3[1], point.shape3[1+j])
          
          plot.list[[(i-1)*n+j]] = compare.gmwm(obj_list[[i]], obj_list[[j]], split = F, CI = F, plot.emp.wv = F, auto.label.wvar = F,
                                                line.color = line.color, line.type = line.type, point.size = point.size,
                                                point.shape = point.shape, 
                                                title = new.title, title.size = facet.label.size,
                                                axis.label.size = facet.label.size, axis.tick.size = axis.tick.size,
                                                axis.x.label = new.axis.x.label, axis.y.label = new.axis.y.label) +
            commonSet
          
        }
        
        #|1|2|
        #|3|4|
        ### remove axis tick and axis text
        if(i == 1 && j ==1){
          #1
          plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + keep.title + remove.x + keep.y
          
        }else if(i == 1 && j == n){
          #2
          plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + keep.title + remove.x + remove.y
          
        }else if(i == n && j == 1){
          #3
          plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + keep.x + keep.y
          
        }else if(i == n && j == n){
          #4
          plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + keep.x + remove.y
        }
        
        #| |--1--| |
        #-----------
        #| |     | |
        #|2|--3--|4|
        #| |     | |
        #-----------
        #| |--5--| |
        
        if(n > 2){
          if (i == 1 && j %in% 2:(n-1)){
            #1
            plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + keep.title + remove.x + remove.y
            
          }else if ( i%in% 2:(n-1) && j == 1){
            #2
            plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + remove.x + keep.y
            
          }else if (i %in% 2:(n-1) && j %in% 2:(n-1) ){
            #3
            plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + remove.x + remove.y
            
          }else if(i %in% 2:(n-1) && j == n){
            #4
            plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + remove.x + remove.y
            
          }else if (i == n && j%in% 2:(n-1)){
            #5
            plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] + remove.title + keep.x + remove.y
            
          }
        }
        
        plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] +
          theme(plot.margin = unit(c(0,0,-0.8,0), "cm"))
        
        #         #remove the space between each graph
        #         if(i==1 && j==1){
        #           plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] +
        #             theme(plot.margin = unit(c(0,0,0,0), "cm") )
        #         }else if(i==1 && j%in%2:n ){
        #          plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] +
        #            theme(plot.margin = unit(c(0,0,0,-0.5), "cm") )
        #         }else if(i%in%2:n && j==1){
        #           plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] +
        #             theme(plot.margin = unit(c(-0.8,0,0,0), "cm") )
        #         }else{
        #           plot.list[[(i-1)*n+j]] = plot.list[[(i-1)*n+j]] +
        #             theme(plot.margin = unit(c(-0.8,0,0,-0.5), "cm") )
        #         }
        
        if(i==1&&j==1){
          y.range.low = ggplot_build(plot.list[[(i-1)*n+j]])$panel$ranges[[1]]$y.range[1]
          y.range.high = ggplot_build(plot.list[[(i-1)*n+j]])$panel$ranges[[1]]$y.range[2]
        }else{
          y.range.low = min(y.range.low, ggplot_build(plot.list[[(i-1)*n+j]])$panel$ranges[[1]]$y.range[1])
          y.range.high = max(y.range.high, ggplot_build(plot.list[[(i-1)*n+j]])$panel$ranges[[1]]$y.range[2])
        }
        
      }#end of j loop
    }#end of i loop
    
    for(i in 1:n^2){
      plot.list[[i]] = plot.list[[i]] + coord_cartesian(ylim = c(10^(y.range.low), 10^(y.range.high))) 
      #plot.list[[i]] = plot.list[[i]] + coord_cartesian(ylim = c(1E-10, 1E-4) )
    }
    
    #grid.arrange is just one wrapper function of arrangeGrob
    args.list <- c(plot.list,list(nrow=n,ncol=n, left = textGrob(label = axis.y.label, rot = 90, gp=gpar(fontsize=axis.label.size) ),
                                  right = textGrob(label = ' ', rot=90),
                                  top = textGrob(title, gp=gpar(fontsize=title.size)),
                                  bottom = textGrob(label = axis.x.label, vjust = 0.8, gp=gpar(fontsize=axis.label.size))) )
    do.call("grid.arrange", args.list)
  }
}