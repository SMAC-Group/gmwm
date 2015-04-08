#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param data A \code{matrix} or \code{data.frame} object with only column (e.g. \eqn{N \times 1}{ N x 1 })
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU or SSM
#' @param compute.v A \code{string} indicating the type of covariance matrix solver. "fast", "bootstrap", "asymp.diag", "asymp.comp", "fft"
#' @param augmented A \code{boolean} indicating whether to add additional moments (e.g. mean for drift and variance for all other components).
#' @param p A \code{double} between 0 and 1 that correspondings to the \eqn{\frac{\alpha}{2}}{alpha/2} value for the wavelet confidence intervals.
#' @param robust A \code{boolean} indicating whether to use the robust computation (TRUE) or not (FALSE).
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency.
#' @param G An \code{integer} to sample the space for IMU and SSM models to ensure optimal identitability.
#' @param K An \code{integer} that controls how many times the bootstrapping procedure will be initiated.
#' @param H An \code{integer} that indicates how many different samples the bootstrap will be collect.
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
#' 2. Guessing ARMA models is NOT supported. You must supply specific parameters (see example below)
#' 3. ARMA Bootstrapping.
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
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' data = gen.ts(AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 1), n)
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
gmwm = function(model, data, model.type="imu", compute.v="fast", inference = FALSE, augmented=FALSE, p = 0.025, robust=FALSE, eff=0.6, G=1000, K=1, H = 100){
  
  # Are we receiving one column of data?
  if( (class(data) == "data.frame" && ncol(data) > 1) || ( class(data) == "matrix" && ncol(data) > 1 ) ){
    stop("The function can only operate on one column of data")
  }
  
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
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
  

  # Needed if model contains a drift. 
  
  theta = model$theta

  out = .Call('GMWM_gmwm_master_cpp', PACKAGE = 'GMWM', data, theta, desc, obj, model.type, starting = model$starting,
                                                         p = p, compute_v = compute.v, K = K, H = H, G = G,
                                                         robust=robust, eff = eff, inference = inference)
  
  #colnames(out) = model$desc
  
  estimate = out[[1]]
  rownames(estimate) = model$process.desc
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  
  out = structure(list(estimate = estimate,
                       init.guess = init.guess,
                       wv.empir = out[[3]], 
                       ci.low = out[[4]], 
                       ci.high = out[[5]],
                       compute.v = compute.v,
                       augmented = augmented,
                       V = out[[6]],
                       orgV = out[[7]],                       
                       scales = scales, 
                       robust = robust,
                       eff = eff,
                       model.type = model.type,
                       expect.diff = out[[8]],
                       N = N,
                       G = G,
                       H = H,
                       K = K,
                       theo = out[[9]],
                       decomp.theo = out[[10]],
                       model = model,
                       starting = model$starting,
                       inference = inference,
                       ci.theta = out[[11]],
                       gof.test = out[[12]]), class = "gmwm")
  invisible(out)
}

#' @title Update GMWM object for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param object A \code{gmwm} object.
#' @param model A \code{ts.model} object containing one of the allowed models.
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
#' data = gen.ts(exact.model)
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
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(np > length(object$scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
  }
  
  if(object$robust){
    np = np+1
    if(np > length(object$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need one additional scale since robust requires the amount of parameters + 1 to estimate.")
    }
  }
  
  out = .Call('GMWM_gmwm_update_cpp', PACKAGE = 'GMWM',
                  model$theta,
                  desc, obj, 
                  object$model.type, object$N, object$expect.diff, 
                  object$orgV, object$scales, object$wv.empir,
                  model$starting, 
                  object$compute.v, object$K, object$H,
                  object$G, 
                  object$robust, object$eff)

  estimate = out[[1]]
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

#' @title Summary of GMWM object
#' @description Displays summary information about GMWM object
#' @method summary gmwm
#' @param x A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen.ts(AR1(phi=.1, sigma2 = 1) + AR2(phi=0.95, sigma2 = .1),n)
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' summary(mod)
summary.gmwm = function(object, ...){
  cat("Parameter estimates: \n")
  out = as.matrix(object$estimate)
  colnames(out) = "Estimates"
  print(out)
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
#' mod = gmwm(AR1(), data=x, model.type="imu")
#' plot(mod)
plot.gmwm = function(x, individual = FALSE, ...){
  if (individual == TRUE){
    class(x) = "gmwm2"
    autoplot(x)
  }else{
    autoplot(x)
  }
}

#' @title Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM for each latent process.
#' @method autoplot gmwm2
#' @param object A \code{GMWM2} object.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM for each latent process.
#' @author JJB
autoplot.gmwm2 = function(object, ...){
  .x=low=high=trans_breaks=trans_format=math_format=NULL
  
  # Find number of latent processes
  L = length(object$model$desc) + 1
  
  # Get names of latent processes
  nom = letters[1:L]
  
  # Construct data.frame
  df = data.frame(scales = rep(object$scales,L), WV = c(as.vector(object$decomp.theo),apply(object$decomp.theo,1,sum)), process = rep(nom, each = length(object$scales)))
  WV = data.frame(var = object$wv.empir, low = object$ci.low, high = object$ci.high, scale = object$scales)
    
  CI = ggplot(WV, aes(x = scale, y = low), colour = "#003C7D") + geom_line(linetype = "dotted", colour = "#003C7D") +
    geom_line(aes(y = high),linetype = "dotted", colour = "#003C7D") +
    geom_line(aes(y = var)) + geom_point(aes(y = var), size = 3) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1, fill = "#003C7D") +
    geom_line(aes(x = df$scales, y = df$WV, color = rep(nom, each = length(object$scales))))  +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
          legend.justification=c(0,0), legend.position=c(0,0)) +
    scale_colour_discrete(name  =" ", labels=c(object$model$desc,paste(object$model$desc, collapse = ' + '))) +
    coord_cartesian(ylim = c(min(object$ci.low), (1.05)*max(object$ci.high)))
  
  CI
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
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' mod = gmwm(2*AR1(), data = x)
#' autoplot(mod)
autoplot.gmwm = function(object, ...){
  
  low=high=emp=theo=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
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
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
autoplot.comp = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#EF8A1C","LINE4"="#6E8BBF")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
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
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = TRUE)}
autoplot.compSplit = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
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
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
compare.models = function(GMWM1, GMWM2, split = FALSE){
  x = data.frame(wv.empir = GMWM1$wv.empir, ci_low = GMWM1$ci.low, 
                 ci_high = GMWM1$ci.high, scales = GMWM1$scales, theo1 = GMWM1$theo, theo2 = GMWM2$theo) 
  if (split == TRUE){
    class(x) = "compSplit"
  }else{
    class(x) = "comp"
  }
  autoplot(x)
}
