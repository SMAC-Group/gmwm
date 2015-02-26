

#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @usage gmwm(model, wvcov, signal, model.type="imu", params = list(ar = 1, ma = 2), B = 1000)
#' @param model A \code{ts.model} object containing one of the allowed models. Not used if \code{type="arma"}
#' @param wvcov A \code{wvcov} object
#' @param signal A \code{vec} that is the time series being studied
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU, ARMA, or SSM
#' @param params A \code{list} containing the numbers of AR and MA parameters. Only used if \code{type="arma"}
#' @param B A \code{integer} to sample the space for IMU and SSM models to ensure AR1 identitability.
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
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' save2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#'  
#' # ARMA case
#' set.seed(1336)
#' x = arima.sim(n = 200, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796))
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save = gmwm(wvcov=out, signal=x, model.type="arma", params=list(ar=2,ma=2))
gmwm = function(model, wvcov, signal, model.type="imu", params = list(ar = 1, ma = 2), B = 1000){
  
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1, DR, RW, QN, WN). Do NOT use params!")
  }
  desc = model$desc
  
  if(model.type == "imu" || model.type == "ssm"){
    if(length(desc) == 0 ){
      stop("desc must contain a list of components (e.g. AR1, DR, RW, QN, WN). Do NOT use params!")
    }
    
    np = num_model_params(desc)
    
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }

    if(wvcov$robust){
      np = np+1
      if(np > length(wvcov$scales)){
        stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
      }
      ind = 1:np
      wvcov$scales = wvcov$scales[ind]
      wvcov$wv.empir = wvcov$wv.empir[ind]
      wvcov$V = wvcov$V[ind,ind]
      wvcov$ci_low = wvcov$ci_low[ind]
      wvcov$ci_high = wvcov$ci_high[ind]
    }
    
    N = length(signal)
    
    out = .Call('GMWM_gmwm_imu_ssm_cpp', PACKAGE = 'GMWM', desc, signal, model.type, wvcov$V, wvcov$wv.empir, wvcov$scales, N, B)

    theo = theo_wv(out, desc, wvcov$wv.empir, wvcov$scales, N)
    
  } else if(model.type == "arma"){
    if(length(params) != 2){
      stop("params must be a list with a length of two. The AR component is first and the MA component is second.")
    }
    p = params[[1]]
    q = params[[2]]
    np = p+q+1
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }

    if(wvcov$robust){
      np = np+1
      if(np > length(wvcov$scales)){
        stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
      }
      ind = 1:np
      wvcov$scales = wvcov$scales[ind]
      wvcov$wv.empir = wvcov$wv.empir[ind]
      wvcov$V = wvcov$V[ind,ind]
      wvcov$ci_low = wvcov$ci_low[ind]
      wvcov$ci_high = wvcov$ci_high[ind]
    }
    
    guess.start = arima(signal, order=c(p,0,q), include.mean=FALSE, method="CSS")
    theta = c(coef(guess.start), guess.start$sigma2)
    
    out = .Call('GMWM_gmwm_arma_cpp', PACKAGE = 'GMWM', theta, wvcov$V, p,  q, wvcov$scales, wvcov$wv.empir)
    

    if(p == 0){
      ar = numeric()
    }else{
      ar = out[1:p]
    }
    
    if(q == 0){
      ma = numeric()
    }else{
      ma = out[(p+1):(length(out)-1)]
    }

    theo = arma_to_wv(ar, ma, wvcov$scales, out[length(out)])
  } 
    
  out = structure(list(estimate = as.vector(out),
                       wv.empir=wvcov$wv.empir, 
                       ci_low=wvcov$ci_low, 
                       ci_high = wvcov$ci_high, 
                       scales = wvcov$scales, 
                       theo = theo, 
                       robust = wvcov$robust,
                       eff = wvcov$eff,
                       type = model.type), class = "gmwm")
  invisible(out)
}






#
#adv.gmwm = function(desc, theta, wvcov){
#  .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, wvcovdesc, V, wv_empir, tau, N)
#}


#' @title Wrapper to Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method plot gmwm
#' @param object A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' mod = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
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
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' mod = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
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
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = FALSE)
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
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = TRUE)
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
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = FALSE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' GMWM1 = gmwm(AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' GMWM2 = gmwm(2*AR1(), wvcov=out, signal=x, model.type="imu", B = 10000)
#' compare.models(GMWM1, GMWM2, split = FALSE)
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
